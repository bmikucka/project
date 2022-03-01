#!/usr/bin/perl -s
#*************************************************************************
#
#   Program:    xvalidate
#   File:       xvalidate.pl
#   
#   Version:    V1.0
#   Date:       24.03.22
#   Function:   
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
#   This program is not in the public domain, but it may be copied
#   according to the conditions laid out in the accompanying file
#   COPYING.DOC
#
#   The code may be modified as required, but any modifications must be
#   documented so that the person responsible can be identified. If 
#   someone else breaks this code, I don't want to be blamed for code 
#   that does not work! 
#
#   The code may not be sold commercially or included as part of a 
#   commercial product except as described in the file COPYING.DOC.
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
#   V1.0   24.02.22  Original   By: ACRM
#
#*************************************************************************
# Add the path of the executable to the library path
use FindBin;
use lib $FindBin::Bin;
# Or if we have a bin directory and a lib directory
#use Cwd qw(abs_path);
#use FindBin;
#use lib abs_path("$FindBin::Bin/../lib");

UsageDie() if(defined($::h));

# Default to 10x x-validation
$::n = 10 if(!defined($::n));

# Read in the CSV data file
my $ndata  = 0;
my $header = <>;
while(<>)
{
    push @data, $_;
    $ndata++;
}

my @boundaries = CalculateBoundaries($ndata, $::n);
ShiftBoundaries(\@boundaries, $::n-1, \@data, $ndata);
WriteFolds(\@boundaries, $::n-1, \@data, $ndata, $header);

#*************************************************************************
sub WriteFolds
{
    my($aBoundaries, $nBoundaries, $aData, $ndata, $header) = @_;
    my $lowerBound = 0;
    my $i;
    for($i=0; $i<$nBoundaries; $i++)
    {
        my $fnm = sprintf("fold_%02d.csv", $i);
        if(open(my $fp, '>', $fnm))
        {
            print $fp $header;
            for(my $record=$lowerBound;
                $record < $$aBoundaries[$i];
                $record++)
            {
                print $fp $$aData[$record];
            }
            $lowerBound = $$aBoundaries[$i];
            close($fp);
        }
    }
    my $fnm = sprintf("fold_%02d.csv", $i);
    if(open(my $fp, '>', $fnm))
    {
        print $fp $header;
        for(my $record=$lowerBound;
            $record < $ndata;
            $record++)
        {
            print $fp $$aData[$record];
        }
        close($fp);
    }
}    

#*************************************************************************
sub ShiftBoundaries
{
    my($aBoundaries, $nBoundaries, $aData, $ndata) = @_;

    for(my $i=0; $i<$nBoundaries; $i++)
    {
        # Initial boundary
        my $boundary = $$aBoundaries[$i];

        # Look ahead in the list to find the UP code change
        my $boundaryUp = 0;
        while(UniProtMatch($$aData[$boundary+$boundaryUp],
                           $$aData[$boundary+$boundaryUp+1]))
        {
            $boundaryUp++;
        }

        # Look back in the list to find the UP code change
        my $boundaryDown = 0;
        while(UniProtMatch($$aData[$boundary-$boundaryDown],
                           $$aData[$boundary-$boundaryDown-1]))
        {
            $boundaryDown++;
        }

        # Update the boundary to the smaller shift
        if($boundaryUp < $boundaryDown)
        {
            $$aBoundaries[$i] = $boundary + $boundaryUp + 1;
        }
        else
        {
            $$aBoundaries[$i] = $boundary - $boundaryDown;
        }
    }
}

#*************************************************************************
sub UniProtMatch
{
    my($datum1, $datum2) = @_;
    my @fields1 = split(/:/, $datum1);
    my @fields2 = split(/:/, $datum2);

    if($fields1[1] eq $fields2[1])
    {
        return(1);
    }
    return(0);
}



#*************************************************************************
sub CalculateBoundaries
{
    my($ndata, $nfolds) = @_;
    
    # Calculate chunk size
    my $chunk = int($ndata/$nfolds);

    # Calculate the boundaries
    my @boundaries = ();
    for(my $i=0; $i<($nfolds-1); $i++)
    {
        $boundaries[$i] = $chunk * ($i + 1);
    }
    return(@boundaries);
}

#*************************************************************************
sub UsageDie
{
    print <<__EOF;

xvalidate V1.0 (c) 2022 UCL, Prof. Andrew C.R. Martin

Usage: xvalidate [-x=n] in.csv
       -x Specify number of cross-validation folds [Default: 10]

Takes a .csv file in the format used by our json2csv scripts from the
SAAP pipeline and splits into n separate files ensuring that each UniProt
accession only falls into one file.
    
__EOF
    exit 0;
}
