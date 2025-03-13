#!/usr/bin/perl
use strict;
use warnings;

# Define the necessary paths and URLs
my $base_path = '/mnt/storagelv/home/inesbrancomartins/Tese/tool1_fast';
my $taxonomy_files_dir = "$base_path/taxonomy_files";
my $scripts_dir = "$base_path/scripts";

# Create directories if they do not exist
mkdir $taxonomy_files_dir unless -d $taxonomy_files_dir;
mkdir $scripts_dir unless -d $scripts_dir;

# URL for taxonomy files
my $taxdmp_url = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip";

# Download taxonomy files
print "Downloading taxonomy files...\n";
system("wget -q -O $taxonomy_files_dir/taxdmp.zip $taxdmp_url");

# Unzip the downloaded file
print "Unzipping taxonomy files...\n";
system("unzip -q $taxonomy_files_dir/taxdmp.zip -d $taxonomy_files_dir");

# Execute the Python script
print "Executing Python script...\n";
system("python3 $scripts_dir/taxonomy_hierarchy.py");

print "Configuration completed.\n";
