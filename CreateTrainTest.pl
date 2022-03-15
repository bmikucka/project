#!/usr/bin/perl

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
