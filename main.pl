#!/usr/bin/perl
use strict;
use warnings;
use Time::HiRes qw(time);
use FindBin qw($Bin);
use File::Spec::Functions qw(catdir catfile);

# User-modifiable parameters
my $mash_threshold_refseq = 0.90;  # Initial threshold for RefSeq
my $mash_threshold_gtdb = 0.90;    # Initial threshold for GTDB
my $mash_threshold_custom = 0.90;  # Initial threshold for the custom database

my $classification_processes = 8;  # Default value: 4

# Define the base path as the directory containing this script
my $base_path = $Bin;

# Define paths for the scripts
my $scripts_dir = catdir($base_path, 'scripts');
my $mash_script = catfile($scripts_dir, 'mash.sh');
my $download_script = catfile($scripts_dir, 'downloadDB.py');
my $minimap_script = catfile($scripts_dir, 'minimap2.sh');
my $classification_script = catfile($scripts_dir, 'classification.py');
my $cleandf_script = catfile($scripts_dir, 'cleandf.py');

# Define common paths
my $output_dir = catdir($base_path, 'output');
my $data_dir = catdir($base_path, 'data');
my $cache_dir = catdir($base_path, 'cache');

# Prompt the user for the input directory (where the .fna files are located)
print "Please enter the path to the input directory (containing .fna files): ";
chomp(my $input_dir = <STDIN>);

# Validate the input directory
unless (-d $input_dir) {
    die "Error: The input directory '$input_dir' does not exist.\n";
}

# Specific paths for each script
my %paths = (
    mash => {
        mash_screen => catfile($data_dir, 'sketch1.msh'),
        gtdb_mash_screen => catfile($data_dir, 'sketch2.msh'),
        custom_mash_screen => catfile($data_dir, 'sketch3.msh'),  # Path for the third database
        screen_tab => catfile($output_dir, 'screen.tab'),
        gtdb_screen_tab => catfile($output_dir, 'gtdb_screen.tab'),
        custom_screen_tab => catfile($output_dir, 'custom_screen.tab'),  # Output file for the third database
        filtered_screen => catfile($output_dir, 'filtered_screen.tab'),
        sorted_screen => catfile($output_dir, 'sorted_screen.tab'),
        top_hits => catfile($output_dir, 'top_hits.tab'),
        selected_genomes => catfile($output_dir, 'selected_genomes.txt'),
    },
    download => {
        genomes_file => catfile($output_dir, 'selected_genomes.txt'),
        downloaded_genomes => catdir($data_dir, 'downloaded_genomes'),
        taxonomy_file => catfile($data_dir, 'detailed_taxonomy.tsv'),
        cache_dir => $cache_dir,
        log_file => catfile($data_dir, 'downloaded_genomes', 'genome_download_log.txt'),
    },
    minimap => {
        reference_set => catfile($data_dir, 'downloaded_genomes', 'combined_genomes.fasta'),
        nt_mmi => catfile($output_dir, 'reference.mmi'),
        resultados_paf => catfile($output_dir, 'resultados.paf'),
    },
    classification => {
        paf_file => catfile($output_dir, 'resultados.paf'),
        taxonomy_file => catfile($data_dir, 'detailed_taxonomy.tsv'),
        hierarchy_file => catfile($data_dir, 'taxonomy_hierarchy.tsv'),
        output_file => catfile($output_dir, 'classified_sequences.tsv'),
    },
);

# Function to run a command and handle errors
sub run_command {
    my ($command) = @_;
    print "Executing: $command\n";
    my $output = `$command 2>&1`;
    my $exit_status = $? >> 8;

    if ($exit_status == 0) {
        print "Success: $output\n";
    } else {
        print "Error executing: $output\n";
        exit($exit_status);
    }
}

# Start timing the execution
my $start_time = time();

# Step 1: Mash with RefSeq dataset
run_command("$mash_script $input_dir $paths{mash}{mash_screen} $paths{mash}{screen_tab} $paths{mash}{filtered_screen} $paths{mash}{sorted_screen} $paths{mash}{top_hits} $paths{mash}{selected_genomes} $mash_threshold_refseq");

# Step 2: Mash with GTDB dataset
run_command("$mash_script $input_dir $paths{mash}{gtdb_mash_screen} $paths{mash}{gtdb_screen_tab} /tmp/gtdb_filtered.tab /tmp/gtdb_sorted.tab /tmp/gtdb_top_hits.tab /tmp/gtdb_selected_genomes.txt $mash_threshold_gtdb");
run_command("cat /tmp/gtdb_selected_genomes.txt >> $paths{mash}{selected_genomes}");
run_command("sort -u -o $paths{mash}{selected_genomes} $paths{mash}{selected_genomes}");

# Step 3: Mash with Custom dataset
run_command("$mash_script $input_dir $paths{mash}{custom_mash_screen} $paths{mash}{custom_screen_tab} /tmp/custom_filtered.tab /tmp/custom_sorted.tab /tmp/custom_top_hits.tab /tmp/custom_selected_genomes.txt $mash_threshold_custom");
run_command("cat /tmp/custom_selected_genomes.txt >> $paths{mash}{selected_genomes}");
run_command("sort -u -o $paths{mash}{selected_genomes} $paths{mash}{selected_genomes}");

# Step 4: Download genomes
run_command("python3 $download_script $paths{download}{genomes_file} $paths{download}{downloaded_genomes} $paths{download}{taxonomy_file} $paths{download}{cache_dir}");

# Step 5: Run Minimap2
run_command("$minimap_script $input_dir $paths{minimap}{reference_set} $paths{minimap}{nt_mmi} $paths{minimap}{resultados_paf}");

# Step 6: Classification
run_command("python3 $classification_script --paf $paths{classification}{paf_file} --taxonomy $paths{classification}{taxonomy_file} --hierarchy $paths{classification}{hierarchy_file} --output $paths{classification}{output_file} --processes $classification_processes");

# Calculate and display execution time
my $end_time = time();
my $execution_time = $end_time - $start_time;

print "All steps were executed successfully.\n";
print "\nExecution Summary:\n";
print "Total execution time: " . sprintf("%.2f", $execution_time) . " seconds\n";