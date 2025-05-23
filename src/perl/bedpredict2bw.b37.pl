#!/usr/bin/perl
use strict;
use warnings;

# Alert user and halt if error occurs while the perl script is executed
sub run_command {
    my ($cmd) = @_;
	# Capture STDOUT and STDERR
    my $captured_output = `$cmd 2>&1`;
    my $exit_code = $? >> 8;
    my $signal_num = $? & 127;
    my $dumped_core = $? & 128;

    if ($exit_code != 0 || $signal_num != 0 || $dumped_core) {
        my $error_message = "Error: Command execution failed.\n";
        $error_message .= "  Command: $cmd\n";
        $error_message .= "  Exit code: $exit_code\n";
        $error_message .= "  Signal: $signal_num\n" if $signal_num;
        $error_message .= "  Core dumped\n" if $dumped_core;
        $error_message .= "  Output (STDERR/STDOUT):\n$captured_output\n" if $captured_output;
        die $error_message;
    }
    return $captured_output; 
}


if (@ARGV != 3) {
    die "Usage: $0 <output_prefix> <input_file> <genome_sizes_file>\n" .
        "  <output_prefix>:     Prefix for all output files (e.g., 'my_sample').\n" .
        "  <input_file>:        Path to the input methylation data file (e.g., methylation_calls.txt or .txt.gz).\n" .
        "  <genome_sizes_file>: Path to the chromosome sizes file for the reference genome (e.g., hg19.chrom.sizes).\n";
}
my ($output_prefix, $input_file, $genome_sizes_file) = @ARGV;

print "INFO: Starting processing for output prefix '$output_prefix'.\n";
print "INFO: Input file: '$input_file'.\n";
print "INFO: Genome sizes file: '$genome_sizes_file'.\n";


die "Error: Input file '$input_file' not found.\n" unless -e $input_file;
die "Error: Genome sizes file '$genome_sizes_file' not found.\n" unless -e $genome_sizes_file;


my $methy_bedgraph = "$output_prefix.methy.bedgraph";
my $cov_bedgraph = "$output_prefix.cov.bedgraph";
my $m_count_bedgraph = "$output_prefix.methy_count.bedgraph";


my $methy_bw = "$output_prefix.methy.b37.bw"; 
my $cov_bw = "$output_prefix.cov.b37.bw";
my $m_count_bw = "$output_prefix.methy_count.b37.bw";

print "INFO: Generating intermediate BedGraph files...\n";

my $initial_processor_cmd;
my $post_processor_pipe_cmd = "sed 's/^chrM/MT/g;s/^chr//g;'";

if ($input_file =~ /\.gz$/) {
    print "INFO: Input file is gzipped. Using zcat.\n";
    $initial_processor_cmd = "zcat '$input_file' | grep -v \"start\"";
} else {
    print "INFO: Input file is not gzipped. Using grep on file.\n";
    $initial_processor_cmd = "grep -v \"start\" '$input_file'";
}

my @bedgraph_cmds_to_run;
push @bedgraph_cmds_to_run, "$initial_processor_cmd | $post_processor_pipe_cmd | cut -f1-3,5 | bedtools sort -i stdin > '$m_count_bedgraph'";
push @bedgraph_cmds_to_run, "$initial_processor_cmd | $post_processor_pipe_cmd | cut -f1-3,6 | bedtools sort -i stdin > '$cov_bedgraph'";
push @bedgraph_cmds_to_run, "$initial_processor_cmd | $post_processor_pipe_cmd | cut -f1-4   | bedtools sort -i stdin > '$methy_bedgraph'";

foreach my $cmd (@bedgraph_cmds_to_run) {
    run_command($cmd);
}

print "INFO: Intermediate BedGraph files generated.\n";

# Convert BedGraph to BigWig
print "INFO: Converting BedGraph files to BigWig format...\n";

die "Error: Methy BedGraph file '$methy_bedgraph' was not created or is empty. Cannot generate BigWig.\n" 
    unless -s $methy_bedgraph; 
run_command("bedGraphToBigWig '$methy_bedgraph' '$genome_sizes_file' '$methy_bw'");
print "INFO: Created BigWig: $methy_bw\n";

die "Error: Coverage BedGraph file '$cov_bedgraph' was not created or is empty. Cannot generate BigWig.\n" 
    unless -s $cov_bedgraph;
run_command("bedGraphToBigWig '$cov_bedgraph' '$genome_sizes_file' '$cov_bw'");
print "INFO: Created BigWig: $cov_bw\n";

die "Error: Methylated Count BedGraph file '$m_count_bedgraph' was not created or is empty. Cannot generate BigWig.\n" 
    unless -s $m_count_bedgraph;
run_command("bedGraphToBigWig '$m_count_bedgraph' '$genome_sizes_file' '$m_count_bw'");
print "INFO: Created BigWig: $m_count_bw\n";

print "INFO: BigWig files generated.\n";

print "INFO: Cleaning up intermediate BedGraph files...\n";
my @files_to_unlink = ($methy_bedgraph, $cov_bedgraph, $m_count_bedgraph);
my $unlinked_count = 0;
foreach my $file_to_unlink (@files_to_unlink) {
    if (-e $file_to_unlink) {
        if (unlink $file_to_unlink) {
            $unlinked_count++;
        } else {
            warn "Warning: Could not unlink intermediate file '$file_to_unlink': $!\n";
        }
    }
}

if ($unlinked_count > 0) {
    print "INFO: $unlinked_count intermediate file(s) unlinked.\n";
} else {
    print "INFO: No intermediate files needed unlinking or were found to unlink.\n";
}

print "INFO: Script finished successfully.\n";