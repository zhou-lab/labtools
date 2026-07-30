import os
import pysam

def split_bam_to_chunks(bam_file, chunk_size, output_dir):
    # Ensure output folder exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Get the base name of the original BAM file
    bam_base_name = os.path.basename(bam_file).replace('.bam', '')
    
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        output_bam = None
        current_chunk_size = 0
        chunk_number = 1
        
        for read in bam:
            if current_chunk_size == 0:
                # Start a new chunk
                output_bam_path = os.path.join(output_dir, f"{bam_base_name}_tmp{chunk_number}.bam")
                output_bam = pysam.AlignmentFile(output_bam_path, "wb", template=bam)
            
            output_bam.write(read)
            current_chunk_size += 1
            
            if current_chunk_size >= chunk_size:
                output_bam.close()
                current_chunk_size = 0
                chunk_number += 1
        
        # Close the last BAM file if it's still open
        if output_bam is not None:
            output_bam.close()