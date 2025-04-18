import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import os

def phred_to_score(char):
    """Convert Phred character to quality score"""
    return ord(char) - 33

def score_to_phred(score):
    """Convert quality score to Phred character"""
    return chr(score + 33)

def trim_sequence_base(record, base_threshold):
    """Trim low-quality bases from both ends (base-by-base approach)"""
    quals = record.letter_annotations["phred_quality"]
    
    # Find left trim position
    left = 0
    while left < len(quals) and quals[left] < base_threshold:
        left += 1
    
    # Find right trim position
    right = len(quals) - 1
    while right >= 0 and quals[right] < base_threshold:
        right -= 1
    
    if left > right:  # Entire sequence is low quality
        return 0, 0, None
    
    # Create new SeqRecord with trimmed sequence and qualities
    trimmed_seq = record.seq[left:right+1]
    trimmed_qual = quals[left:right+1]
    
    trimmed_record = SeqRecord(
        Seq(str(trimmed_seq)),
        id=record.id,
        description=record.description,
        letter_annotations={"phred_quality": trimmed_qual}
    )
    
    return left, len(record) - right - 1, trimmed_record

def trim_sequence_window(record, window_threshold, window_size):
    """Trim low-quality bases using sliding window approach"""
    quals = record.letter_annotations["phred_quality"]
    seq_len = len(quals)
    
    # Calculate window averages from left
    left = 0
    while left <= seq_len - window_size:
        window = quals[left:left+window_size]
        if sum(window) / window_size >= window_threshold:
            break
        left += 1
    
    # Calculate window averages from right
    right = seq_len
    while right >= window_size:
        window = quals[right-window_size:right]
        if sum(window) / window_size >= window_threshold:
            break
        right -= 1
    
    if left > right - window_size:  # No good window found
        return 0, 0, None
    
    # Create new SeqRecord with trimmed sequence and qualities
    trimmed_seq = record.seq[left:right]
    trimmed_qual = quals[left:right]
    
    trimmed_record = SeqRecord(
        Seq(str(trimmed_seq)),
        id=record.id,
        description=record.description,
        letter_annotations={"phred_quality": trimmed_qual}
    )
    
    return left, seq_len - right, trimmed_record

def process_fastq(input_file, base_threshold=None, window_threshold=None, window_size=5):
    """Process FASTQ file with the given trimming parameters"""
    stats = defaultdict(int)
    stats['per_sequence'] = []
    output_records = []
    
    # Generate output filename
    base_name = os.path.splitext(input_file)[0]
    output_file = f"{base_name}_trimmed.fastq"
    
    with open(input_file) as handle:
        for record in SeqIO.parse(handle, "fastq"):
            stats['total_sequences'] += 1
            original_length = len(record)
            
            if base_threshold is not None:
                # Apply base-by-base trimming
                left_trim, right_trim, trimmed_record = trim_sequence_base(record, base_threshold)
                method = "base"
            elif window_threshold is not None:
                # Apply window-based trimming
                left_trim, right_trim, trimmed_record = trim_sequence_window(
                    record, window_threshold, window_size
                )
                method = "window"
            
            if trimmed_record is None:
                stats['discarded_sequences'] += 1
                continue
            
            # Update statistics
            trimmed_length = len(trimmed_record)
            stats['total_bases_trimmed'] += (original_length - trimmed_length)
            stats['total_bases_remaining'] += trimmed_length
            
            # Record per-sequence trimming info
            stats['per_sequence'].append({
                'id': trimmed_record.id,
                'original_length': original_length,
                'trimmed_length': trimmed_length,
                'bases_trimmed': original_length - trimmed_length,
                'left_trim': left_trim,
                'right_trim': right_trim,
                'method': method
            })
            
            output_records.append(trimmed_record)
    
    # Write output FASTQ
    with open(output_file, "w") as handle:
        SeqIO.write(output_records, handle, "fastq")
    
    return stats, output_file

def main():
    parser = argparse.ArgumentParser(
        description="Trim low-quality bases from FASTQ files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("input", help="Input FASTQ file")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--base-threshold", type=int, 
                     help="Quality threshold for base-by-base trimming")
    group.add_argument("--window-threshold", type=int, 
                     help="Average quality threshold for window-based trimming")
    parser.add_argument("--window-size", type=int, default=5,
                     help="Window size for window-based trimming (default: %(default)s)")
    
    args = parser.parse_args()
    
    print(f"Processing {args.input} with:")
    if args.base_threshold:
        print(f"  Base-by-base trimming with threshold {args.base_threshold}")
    else:
        print(f"  Window-based trimming with threshold {args.window_threshold} and window size {args.window_size}")
    
    stats, output_file = process_fastq(
        args.input,
        base_threshold=args.base_threshold,
        window_threshold=args.window_threshold,
        window_size=args.window_size
    )
    
    print("\nTrimming statistics:")
    print(f"  Total sequences processed: {stats['total_sequences']}")
    print(f"  Sequences discarded (all low quality): {stats['discarded_sequences']}")
    print(f"  Sequences kept: {stats['total_sequences'] - stats['discarded_sequences']}")
    print(f"  Total bases trimmed: {stats['total_bases_trimmed']}")
    print(f"  Total bases remaining: {stats['total_bases_remaining']}")
    print(f"\nTrimmed sequences saved to {output_file}")
    
    # Print first few sequences' trimming details
    print("\nTrimming details:")
    for seq in stats['per_sequence']:
        print(f"  Sequence {seq['id']}:")
        print(f"    Original length: {seq['original_length']}")
        print(f"    Trimmed length: {seq['trimmed_length']}")
        print(f"    Bases trimmed: {seq['bases_trimmed']}")
        print(f"    Left trim ({seq['method']}): {seq['left_trim']}")
        print(f"    Right trim ({seq['method']}): {seq['right_trim']}")
        print()

if __name__ == "__main__":
    main()