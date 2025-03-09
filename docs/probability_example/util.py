
import os
import subprocess
import pandas as pd
import numpy as np
import ugbio_core.flow_format.flow_based_read as fbr
import random

def preview_sam_file(sam_file, header_lines=20, read_lines=10):
    """
    Prints partial header and the first few reads from a CRAM file.
    """
    # --- 1) Print the first N lines of the header
    cmd_header = f"samtools view -H {sam_file} | head -n {header_lines}"
    print(f"=== HEADER (first {header_lines} lines) ===")
    process = subprocess.Popen(cmd_header, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    header_output, header_err = process.communicate()
    if process.returncode == 0:
        print(header_output.decode('utf-8'))
    else:
        print("Error retrieving header:", header_err.decode('utf-8'))

    # --- 2) Print the first N reads
    cmd_reads = f"samtools view {sam_file} | head -n {read_lines}"
    print(f"=== FIRST {read_lines} READS ===")
    process = subprocess.Popen(cmd_reads, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    reads_output, reads_err = process.communicate()
    if process.returncode == 0:
        print(reads_output.decode('utf-8'))
    else:
        print("Error retrieving reads:", reads_err.decode('utf-8'))


def extract_all_tags_to_dataframe(sam_file, max_reads=None, bp_length=None, aws=False):
    """
    Extracts reads from a CRAM file using samtools view, then parses all 
    standard SAM fields + any optional tags (including multi-value B tags) 
    into a pandas DataFrame.
    """
    env = os.environ.copy()
    
    # (Optional) If you need AWS creds, do that part here. 
    # For brevity, we'll omit it or just copy from your earlier code.

    # Prepare the samtools command
    if max_reads is None:
        command = f"samtools view {sam_file}"
    else:
        command = f"samtools view {sam_file} | head -n {max_reads}"

    # Launch samtools
    process = subprocess.Popen(command, shell=True,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               env=env)
    output, error = process.communicate()

    if process.returncode != 0:
        print("Error running samtools:", error.decode())
        return pd.DataFrame()

    # We'll accumulate rows (dicts) here
    rows = []

    for line in output.decode().splitlines():
        if not line.strip():
            continue
        
        fields = line.split('\t')
        if len(fields) < 11:
            # Malformed or partial line, skip
            continue
        
        # Standard SAM fields
        qname = fields[0]
        flag  = fields[1]
        rname = fields[2]
        pos   = fields[3]
        mapq  = fields[4]
        cigar = fields[5]
        rnext = fields[6]
        pnext = fields[7]
        tlen  = fields[8]
        seq   = fields[9]
        qual  = fields[10]

        # Optionally truncate the read sequence
        if bp_length is not None and seq:
            seq = seq[:bp_length]

        # Create a dict for this read
        row = {
            "QNAME": qname,
            "FLAG":  flag,
            "RNAME": rname,
            "POS":   pos,
            "MAPQ":  mapq,
            "CIGAR": cigar,
            "RNEXT": rnext,
            "PNEXT": pnext,
            "TLEN":  tlen,
            "SEQ":   seq,
            "QUAL":  qual
        }

        # Optional tags start at field[11:]
        for tag_field in fields[11:]:
            # Typically something like: NM:i:0, MD:Z:55, AS:i:55, tp:B:c,0,0,1,1, ...
            parts = tag_field.split(':', 2)
            if len(parts) == 3:
                tag_name, tag_type, tag_val = parts

                # Convert from string using our helper
                parsed_val = parse_tag_value(tag_type, tag_val)

                row[tag_name] = parsed_val
            else:
                # If the tag doesn't conform to "TAG:TYPE:VAL", store raw
                row[tag_field] = None

        rows.append(row)

    # Convert the list of dicts to a DataFrame
    df = pd.DataFrame(rows)
    df['Index'] = np.array([np.int64(x.split('-')[-1]) for x in df['QNAME'].values])
    df['OP'] = np.array([np.int8(int(x.split('_')[1][0])) for x in df['QNAME'].values])

    return df



def Flow_prob_from_sam_record(alignment_file, batch_size=10000, output_file=None, max_hmer_size=20, max_len=600):
    """
    Process an alignment file in batches and save to a single binary file.
    """
    # Initialize output file
    with open(output_file, 'wb') as f:
        batch = []
        for i, record in enumerate(alignment_file):
            # Process each record
            try:
                flow_read = fbr.FlowBasedRead.from_sam_record(record, max_hmer_size=max_hmer_size)
                prob_i = flow_read._flow_matrix.transpose(1, 0)
            except Exception as e:
                print(f"Skipped record {i + 1}: {e}")
                continue

            # Zero-pad or truncate to max_len
            if prob_i.shape[0] > max_len:
                prob_i = prob_i[:max_len,:]
            else:
                padding = ((0, max_len - prob_i.shape[0]), (0, 0))
                prob_i = np.pad(prob_i, padding, mode='constant')


            batch.append(prob_i)

            # If batch is full, write to file and reset
            if (i + 1) % batch_size == 0:
                np_batch = np.array(batch, dtype=np.float32)
                np_batch.tofile(f)
                batch = []

                print(f'Processed {i + 1} records, in batch {i // batch_size}')

        # Write remaining records in the final batch
        if batch:
            np_batch = np.array(batch, dtype=np.float32)
            np_batch.tofile(f)

            print(f'Processed {i + 1} records, in batch {i // batch_size}')

    print(f'Finished processing {i + 1} records')
# 


def parse_tag_value(tag_type, tag_val):
    """
    Convert the tag_val (string) to the appropriate Python object
    depending on tag_type.
    
    For 'B' arrays, parse into a list of values. 
    For 'i' or 'f', parse into int or float. 
    Otherwise return the raw string.
    """
    if tag_type == 'i':
        # Integer
        return int(tag_val)
    elif tag_type == 'f':
        # Float
        return float(tag_val)
    elif tag_type == 'B':
        # B-type array: e.g. 'B:c,0,0,1,1,0' => subtype='c', array='0,0,1,1,0'
        # The first character in tag_val is the subtype (c, C, s, S, i, I, f)
        # Then we have a comma, then the array elements as CSV.
        # Example: 'c,0,0,1,1,0'
        # We'll parse them into a Python list of integers/floats as appropriate.
        #
        # NOTE: The format is: B:<subtype>,<val1>,<val2>,...
        #
        # So let's split on the first comma:
        #   subtype = 'c'
        #   array_part = '0,0,1,1,0'
        
        try:
            subtype, array_part = tag_val.split(',', 1)
        except ValueError:
            # If something is off, just return raw
            return tag_val
        
        # Now parse the array elements:
        elements_str = array_part.split(',')
        
        # Decide how to convert them based on the subtype:
        # c = int8, C = uint8, s = int16, S = uint16, i = int32, I = uint32, f = float
        # For simplicity, we'll treat them as int if subtype != 'f'.
        
        if subtype == 'f':
            # parse as float array
            return [float(x) for x in elements_str]
        else:
            # parse as int array
            return [int(x) for x in elements_str]
    else:
        # By default, return the string as-is (Z, H, A, etc.)
        return tag_val



def key2seq(key, flow_order=None, start=0, return_flow_indices=False):
    """
    Convert a list of counts (key) into a sequence of characters based on the specified flow order.
    
    Parameters
    ----------
    key : list of int
        A list of integers, each indicating how many times a character from the flow order is repeated.
    
    flow_order : list of str, optional
    
    start : int, optional
        The starting offset within the flow order. Default is 0.
    
    return_flow_indices : bool, optional
        If True, the function will return a tuple containing:
          1) The generated sequence (str).
          2) A list of integers (same length as the sequence) where each integer indicates 
             the source flow index of the corresponding character in the sequence.
        Default is False.

    Returns
    -------
    str or tuple
        If return_flow_indices is False (default), returns only the generated sequence (str).
        If return_flow_indices is True, returns a tuple of:
          (generated_sequence_str, list_of_flow_indices).
    """
    if flow_order is None:
        flow_order = ['T', 'G', 'C', 'A']


    seq_chars = []
    flow_indices = []

    for i, count in enumerate(key):
        # Determine the character for this flow position, factoring in the 'start' offset
        current_char = flow_order[(start + i) % len(flow_order)]
        
        # Extend the sequence by 'count' copies of this character
        seq_chars.extend([current_char] * count)
        
        # If needed, track the flow index for each character
        if return_flow_indices:
            flow_index = (start + i) # % len(flow_order)
            flow_indices.extend([flow_index] * count)

    sequence = ''.join(seq_chars)

    if return_flow_indices:
        return sequence, np.array(flow_indices)

    return sequence
