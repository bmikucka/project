#!/usr/bin/perl

#*************************************************************************
#
#   Program:    cleanjson
#   File:       cleanjson.pl
#   
#   Version:    V1.0
#   Function:   Cleans JSON files
#   Author:     Andrew CR Martin
#               
#*************************************************************************
#
#   Description:
#   ============
#
#*************************************************************************
#
#   Usage:
#   ======
#           
#
#*************************************************************************
#
#   Revision History:
#   =================
#   V1.0  Original       By: ACRM
#
#*************************************************************************
use strict;
my $data = '';
my $oldData;

# Read Data
while(<>)
{
    $data .= $_;
}

# Remove repeated commas
do {
    $oldData = $data;
    $data =~ s/,\n,/,/g;
}   while($data ne $oldData);

# Remove trailing commas
do {
    $oldData = $data;
    $data =~ s/,\n(\s+)\]/$1\]/g;
}   while($data ne $oldData);

# Remove leading commas
do {
    $oldData = $data;
    $data =~ s/\[\n,/\[\n/g;
}   while($data ne $oldData);


print $data;
