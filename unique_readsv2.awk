# AWK script to count the number of *distinct* reads that contain at least one occurrence
# of each tandem motif, and calculate the percentage of total unique reads they represent.
# This version skips the first 2 and last 10 lines of the input file.

# Array 'read_names': Stores $1 (Read Name) for all lines. Index is NR.
# Array 'motif_names': Stores $2 (Motif Name) for all lines. Index is NR.
# Array 'count': used for deduping (Motif, Read) pairs.
# Array 'all_reads': used to find the total number of unique reads (the denominator).

{
    # 1. Store the data for all lines in memory.
    # This prepares the data to be processed later, allowing us to reference the total line count (NR).
    read_names[NR] = $3
    motif_names[NR] = $1
}

# The END block runs after all input lines have been processed.
END {
    # Define the range of lines to process: from line 3 up to NR - 10 (inclusive).
    # Lines 1 and 2 are skipped. The last 10 lines (NR-9 through NR) are skipped.
    start_line = 3
    end_line = NR - 10

    # Safety check: ensure there are enough lines to process (i.e., NR >= 13)
    if (end_line < start_line) {
        # Print error message to stderr (>&2) and exit gracefully.
        print "Error: Input file has only " NR " lines. Need at least 13 lines to skip first 2 and last 10." | "cat >&2"
        exit 1
    }

    # Process only the lines within the required range [3, NR - 10]
    for (i = start_line; i <= end_line; i++) {
        read_name = read_names[i]
        motif_name = motif_names[i]

        # Deduplicate by (Motif, Read) pair for per-motif counting
        count[motif_name, read_name] = 1

        # Track all unique read names globally for the denominator calculation
        all_reads[read_name] = 1
    }

    # Determine the total number of unique reads found in the processed range.
    total_unique_reads = length(all_reads)

    # Re-aggregate the distinct read counts per motif from the 'count' array
    for (key in count) {
        # Split the composite key back into its two components: Motif and Read.
        # The key is "MotifSUBSEPRead"
        split(key, parts, SUBSEP)
        final_counts[parts[1]]++
    }

    # Print the final report header
    printf "%-25s %-15s %-10s\n", "TANDEM MOTIF", "DISTINCT READS", "PERCENT (%)"
    print "------------------------- --------------- ----------"

    # Iterate through the final sum array, calculate percentage, and print
    for (motif in final_counts) {
        motif_count = final_counts[motif]

        # Check if total_unique_reads is 0 to avoid division by zero
        if (total_unique_reads > 0) {
            # Calculate percentage: (Motif Reads / Total Reads) * 100
            percentage = (motif_count / total_unique_reads) * 100
        } else {
            percentage = 0.00
        }

        # Print the results, formatting the percentage to two decimal places
        printf "%-25s %15d %10.2f\n", motif, motif_count, percentage
    }

    print "------------------------- --------------- ----------"
    printf "PROCESSED LINE RANGE: %d through %d\n", start_line, end_line
    printf "TOTAL UNIQUE READS (in range): %d\n", total_unique_reads
}
