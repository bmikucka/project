#!/usr/bin/perl

#*************************************************************************
#
#   Program:    CreateTrainTest
#   File:       CreateTrainTest.pl
#   
#   Version:    V1.0
#   Date:       15.03.22
#   Function:   Create test and training subsets for Weka ML
#
#   Copyright:  (c) UCL, Prof. Andrew C. R. Martin, 2022
#   Author:     Prof. Andrew C. R. Martin
#   Address:    Institute of Structural and Molecular Biology
#               Division of Biosciences
#               University College
#               Gower Street
#               London
#               WC1E 6BT
#   EMail:      andrew@bioinf.org.uk
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
#*************************************************************************
#
#   Revision History:
#   =================
#   V1.0   15.03.22  Original   By: ACRM
#
#*************************************************************************


use strict;

my $fold    = shift(@ARGV);
my $nFolds  = shift(@ARGV);
my $foldDir = shift(@ARGV);

# Create the test fold
my $fnmTest = sprintf("test_%02d.csv", $fold);
my $fnmIn   = sprintf("$foldDir/fold_%02d.csv", $fold);
`cp $fnmIn $fnmTest`;

# Create the training fold
my $fnmTrain = sprintf("train_%02d.csv", $fold);
`head -1 $foldDir/fold_00.csv > $fnmTrain`;
for(my $thisFold=0; $thisFold<$nFolds; $thisFold++)
{
    if($thisFold != $fold)
    {
        my $fnmIn  = sprintf("$foldDir/fold_%02d.csv", $thisFold);
        `grep -v Binding $fnmIn >> $fnmTrain`;
    }
}
