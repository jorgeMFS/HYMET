#!/usr/bin/perl
use strict;
use warnings;
use FindBin qw($Bin);
use File::Spec::Functions qw(catdir catfile);

# Use the directory containing this script as the base path
my $base_path = $Bin;

# Define the necessary directories
my $taxonomy_files_dir = catdir($base_path, 'taxonomy_files');
my $scripts_dir = catdir($base_path, 'scripts');
my $data_dir = catdir($base_path, 'data');
my $output_dir = catdir($base_path, 'output');
my $cache_dir = catdir($base_path, 'cache');

# Ensure the data directory exists for downstream scripts
mkdir $data_dir unless -d $data_dir;
mkdir $output_dir unless -d $output_dir;
mkdir $cache_dir unless -d $cache_dir;

# Create directories if they don't exist
mkdir $taxonomy_files_dir unless -d $taxonomy_files_dir;
mkdir $scripts_dir unless -d $scripts_dir;

# URL for taxonomy files
my $taxdmp_url = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip";

# Download taxonomy files
print "Downloading taxonomy files...\n";
my $taxdmp_zip = catfile($taxonomy_files_dir, 'taxdmp.zip');
system('wget', '-q', '-O', $taxdmp_zip, $taxdmp_url);

# Unzip the downloaded file
print "Unzipping taxonomy files...\n";
system('unzip', '-q', $taxdmp_zip, '-d', $taxonomy_files_dir);

# Execute the Python script
print "Executing Python script...\n";
my $taxonomy_script = catfile($scripts_dir, 'taxonomy_hierarchy.py');
system('python3', $taxonomy_script, $taxonomy_files_dir, $data_dir);

print "Configuration completed.\n";