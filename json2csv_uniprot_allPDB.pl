#!/usr/bin/perl
#*************************************************************************
#
#   Program:    json2csv
#   File:       json2csv_unprot_allPDB_pl.txt
#   
#   Version:    V1.0
#   Date:       04.03.22
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
#   V1.0             Original           By: ACRM
#   V1.1   04.03.22  Added SprotFTdist  By: BAM
#
#*************************************************************************
# Import libraries

use strict;
use FindBin;
use Cwd qw(abs_path);
use lib abs_path("$FindBin::Bin/../lib");
use lib abs_path("$FindBin::Bin/");
use JSONSAAP;

my $jsonFile = $ARGV[0];
my $dataset  = $ARGV[1];

my $id= 0;
my ($id,$result, $Binding, $SProtFT, $SprotFTdist, $Interface, $Relaccess, $Impact, $HBonds, $SPhobic, $CPhilic, $BCharge, $SSGeom, $Voids, $MLargest, $NLargest, $Clash, $Glycine, $Proline, $CisPro);

open ( IN, $jsonFile ) || die "Cannot open $jsonFile!\n";

my $content = "";

while(<IN>)
{
    $content .= $_;
}

my $jsonText = JSONSAAP::Decode($content);


my ($type, $error) = JSONSAAP::Check($jsonText);

if($error ne "")
{
    # 06.12.13 Added printing of "ERROR:" to make it easy to check By: ACRM
    print "ERROR: $error\n";
    exit 1;
}

my($uniprotac, $res, $nat, $mut) = JSONSAAP::IdentifyUniprotSaap($jsonText);
my @jsonSaaps = JSONSAAP::GetSaapArray($jsonText);



print  "num:uniprotac:res:nat:mut:pdbcode:chain:resnum:mutation:structuretype:resolution:rfactor,Binding,SProtFT0,SProtFT1,SProtFT2,SProtFT3,SProtFT4,SProtFT5,SProtFT6,SProtFT7,SProtFT8,SProtFT9,SProtFT10,SProtFT11,SProtFT12,SprotFTdist-ACT_SITE,SprotFTdist-BINDING,SprotFTdist-CA_BIND,SprotFTdist-DNA_BIND,SprotFTdist-NP_BIND,SprotFTdist-METAL,SprotFTdist-MOD_RES,SprotFTdist-CARBOHYD,SprotFTdist-MOTIF,SprotFTdist-LIPID,Interface,Relaccess,Impact,HBonds,SPhobic,CPhim,BCharge,SSGeom,Voids,MLargest1,MLargest2,MLargest3,MLargest4,MLargest5,MLargest6,MLargest7,MLargest8,MLargest9,MLargest10,NLargest1,NLargest2,NLargest3,NLargest4,NLargest5,NLargest6,MLargest7,MLargest8,NLargest9,NLargest10,Clash,Glycine,Proline,CisPro,dataset\n";

foreach my $jsonSaaps (@jsonSaaps)
{
    $id++;
    my($file, $pdbcode, $chain, $resnum, $mutation) = JSONSAAP::IdentifyPDBSaap($jsonSaaps);
    my ($structuretype,$resolution,$rfactor) = JSONSAAP::GetPDBExperiment($jsonSaaps);
    print"$id:$uniprotac:$res:$nat:$mut:$pdbcode:$chain:$resnum:$mutation:$structuretype:$resolution:$rfactor,";
    
    my @analyses = JSONSAAP::ListAnalyses($jsonSaaps);
    foreach my $analysis (@analyses)
    {
        my $pResults = JSONSAAP::GetAnalysis($jsonSaaps, $analysis);
        
        if    ($analysis eq "Binding")      { $Binding   = Binding($pResults);}
        elsif ($analysis eq "SProtFT")      { $SProtFT   = SProtFT($pResults);}
        elsif ($analysis eq "SprotFTdist")  { $SprotFTdist = SprotFTdist($pResults);}
        elsif ($analysis eq "Interface")    {($Interface,$Relaccess) = Interface($pResults);}            
        elsif ($analysis eq "Impact")       { $Impact    = Impact($pResults); }
        elsif ($analysis eq "HBonds")       { $HBonds    = HBonds($pResults); }
        elsif ($analysis eq "SurfacePhobic"){ $SPhobic   = SPhobic($pResults);}  
        elsif ($analysis eq "CorePhilic")   { $CPhilic   = CPhilic($pResults);}      
        elsif ($analysis eq "BuriedCharge") { $BCharge   = BCharge($pResults);}    
        elsif ($analysis eq "SSGeom")       { $SSGeom    = SSGeom($pResults); }     
        elsif ($analysis eq "Voids")        {($Voids,$NLargest,$MLargest) = Voids($pResults);}        
        elsif ($analysis eq "Clash")        { $Clash     = Clash($pResults);  }      
        elsif ($analysis eq "Glycine")      { $Glycine   = Glycine($pResults);}     
        elsif ($analysis eq "Proline")      { $Proline   = Proline($pResults);}        
        elsif ($analysis eq "CisPro")       { $CisPro    = CisPro($pResults); }    
    }    
    
    print "$Binding,$SProtFT,$SprotFTdist,$Interface,$Relaccess,$Impact,$HBonds,$SPhobic,$CPhilic,$BCharge,$SSGeom,$Voids,$MLargest,$NLargest,$Clash,$Glycine,$Proline,$CisPro,$dataset\n";    
    
#",\n2\t$Binding,\n3-15\t$SProtFT,\n16\t$Interface,\n17\t$Relaccess,\n18\t$Impact,\n19\t$HBonds,\n20\t$SPhobic,\n21\t$CPhilic,\n22\t$BCharge,\n23\t$SSGeom,\n24\t$Voids,\n25-34\t$NLargest,\n35-44\t$MLargest,\n45\t$Clash,\n46\t$Glycine,\n47\t$Proline,\n48\t$CisPro,\n49\t$dataset\n";
}

#-------------------------------------------------------------------------
sub Binding
{
    my($pResults) = @_;
    my $result;
    
    if ($$pResults{'Binding-BOOL'} ne "")
    {
        if    ($$pResults{'Binding-BOOL'} eq 'OK' ) {$result='0';}
        elsif ($$pResults{'Binding-BOOL'} eq 'BAD') {$result='1';} 
    }
    else 
    {
        $result = '?';
    }   
    return($result);
}        
#-------------------------------------------------------------------------
sub SProtFT
{
    my($pResults) = @_;
    my $result = "";
    
    if ($$pResults{'SProtFT-FEATURES'} ne "")
    {
        $result = $$pResults{'SProtFT-FEATURES'} ;
        $result = substr($result, 0, 1).",".substr($result, 1, 1).",".substr($result, 2, 1).",".substr($result, 3, 1).",".substr($result, 4, 1).",".substr($result, 5, 1).",".substr($result, 6, 1).",".substr($result, 7, 1).",".substr($result, 8, 1).",".substr($result, 9, 1).",".substr($result, 10, 1).",".substr($result, 11, 1).",".substr($result, 12, 1);
    }
    else 
    {
        $result = '?,?,?,?,?,?,?,?,?,?,?,?,?';
    }       
    return($result); 
}
#-------------------------------------------------------------------------
sub SprotFTdist
{
    my($pResults) = @_;

    my $SprotFTdistACT_SITE = '?';
    my $SprotFTdistBINDING  = '?';
    my $SprotFTdistCA_BIND  = '?';
    my $SprotFTdistDNA_BIND = '?';
    my $SprotFTdistNP_BIND  = '?';
    my $SprotFTdistMETAL    = '?';
    my $SprotFTdistMOD_RES  = '?';
    my $SprotFTdistCARBOHYD = '?';
    my $SprotFTdistMOTIF    = '?';
    my $SprotFTdistLIPID    = '?';
    
    if($$pResults{'SprotFTdist-ACT_SITE'} ne "") { $SprotFTdistACT_SITE = $$pResults{'SprotFTdist-ACT_SITE'}; }
    if($$pResults{'SprotFTdist-BINDING'} ne "")  { $SprotFTdistBINDING  = $$pResults{'SprotFTdist-BINDING'};  }
    if($$pResults{'SprotFTdist-CA_BIND'} ne "")  { $SprotFTdistCA_BIND  = $$pResults{'SprotFTdist-CA_BIND'};  }
    if($$pResults{'SprotFTdist-DNA_BIND'} ne "") { $SprotFTdistDNA_BIND = $$pResults{'SprotFTdist-DNA_BIND'}; }
    if($$pResults{'SprotFTdist-NP_BIND'} ne "")  { $SprotFTdistNP_BIND  = $$pResults{'SprotFTdist-NP_BIND'};  }
    if($$pResults{'SprotFTdist-METAL'} ne "")    { $SprotFTdistMETAL    = $$pResults{'SprotFTdist-METAL'};    }
    if($$pResults{'SprotFTdist-MOD_RES'} ne "")  { $SprotFTdistMOD_RES  = $$pResults{'SprotFTdist-MOD_RES'};  }
    if($$pResults{'SprotFTdist-CARBOHYD'} ne "") { $SprotFTdistCARBOHYD = $$pResults{'SprotFTdist-CARBOHYD'}; }
    if($$pResults{'SprotFTdist-MOTIF'} ne "")    { $SprotFTdistMOTIF    = $$pResults{'SprotFTdist-MOTIF'};    }
    if($$pResults{'SprotFTdist-LIPID'} ne "")    { $SprotFTdistLIPID    = $$pResults{'SprotFTdist-LIPID'};    }

    my $result = "$SprotFTdistACT_SITE,$SprotFTdistBINDING,$SprotFTdistCA_BIND,$SprotFTdistDNA_BIND,$SprotFTdistNP_BIND,$SprotFTdistMETAL,$SprotFTdistMOD_RES,$SprotFTdistCARBOHYD,$SprotFTdistMOTIF,$SprotFTdistLIPID";
    
    return($result); 
}

#-------------------------------------------------------------------------
sub Interface
{
    my($pResults) = @_;
    my $result1 = ""; 
    my $result2 = "";
    
    if (($$pResults{'Interface-RELACCESS'} ne "") && ($$pResults{'Interface-RELACCESS-MOL'} ne ""))
    {
        $result1 = $$pResults{'Interface-RELACCESS'} - $$pResults{'Interface-RELACCESS-MOL'};
        # Relative accessibility 
        $result2= $$pResults{'Interface-RELACCESS'};
    }
    else 
    {
        $result1 = '?'; 
        $result2 = '?';
    } 
    return($result1, $result2);
}  
#-------------------------------------------------------------------------
sub Impact
{ 
    my($pResults) = @_;
    my $result = "";
    
    if ($$pResults{'Impact-CONSSCORE'} ne "")
    {
        $result = $$pResults{'Impact-CONSSCORE'} ;
    }
    else 
    {
        $result = '?';
    }
    return($result); 
}
#-------------------------------------------------------------------------
#When HBonds-ENERGY is a number I can't remember if positive or negative is good. Consequently you need to look at the range of numbers and make sure that BAD and OK are assigned appropriate numbers. For example, if a BAD HBond has an energy of (for example) +100 and a good HBond has an energy of -100 then when HBonds-ENERGY is NULL you should use -100 for OK and +100 for BAD.

sub HBonds
{
    my($pResults) = @_;
    my $result = "";
    
    if ($$pResults{'HBonds-ENERGY'} ne "")
    {
        if ($$pResults{'HBonds-ENERGY'} eq "NULL")
        {
            if ($$pResults{'HBonds-BOOL'} eq 'OK')
            {
                $result = '20';
            }
            if ($$pResults{'HBonds-BOOL'} eq 'BAD')
            {
                $result = '-20';
            }
        }
        else 
        {
            $result = $$pResults{'HBonds-ENERGY'};
        }
    }
    else
    {
        $result = '?';
    }        
    
    return($result);
}
#-------------------------------------------------------------------------
sub SPhobic
{
    my($pResults) = @_;           
    my $result = "";
    
    if (($$pResults{'SurfacePhobic-NATIVE-HPHOB'} ne "") && ($$pResults{'SurfacePhobic-MUTANT-HPHOB'} ne ""))
    {
        if (($$pResults{'SurfacePhobic-MUTANT-HPHOB'}) > ($$pResults{'SurfacePhobic-NATIVE-HPHOB'}))
        {
            $result = $$pResults{'SurfacePhobic-NATIVE-HPHOB'} - $$pResults{'SurfacePhobic-MUTANT-HPHOB'}; 
        }
        else
        {
            $result = '0';
        }
    }
    else 
    {
        $result = '?'; 
    }
    return($result);           
}           
#-------------------------------------------------------------------------
sub CPhilic
{
    my($pResults) = @_;
    my $result = "";
    
    if (($$pResults{'CorePhilic-NATIVE-HPHOB'} ne "") && ($$pResults{'CorePhilic-MUTANT-HPHOB'} ne ""))
    {
        if ($$pResults{'CorePhilic-MUTANT-HPHOB'} < $$pResults{'CorePhilic-NATIVE-HPHOB'})
        {
            $result = $$pResults{'CorePhilic-NATIVE-HPHOB'} - $$pResults{'CorePhilic-MUTANT-HPHOB'}; 
        }
        else
        {
            $result = '0'; 
        }    
    }
    else 
    {
        $result = '?';
    }
    return($result);         
}
#-------------------------------------------------------------------------
sub BCharge
{
    my($pResults) = @_;
    my $result="";
    
    if (($$pResults{'BuriedCharge-NATIVE-CHARGE'} ne "") && ($$pResults{'BuriedCharge-MUTANT-CHARGE'} ne ""))
    {
        $result = $$pResults{'BuriedCharge-NATIVE-CHARGE'} - $$pResults{'BuriedCharge-MUTANT-CHARGE'}; 
    }
    else 
    {
        $result = '?';
    }
    return($result);        
}
#-------------------------------------------------------------------------
sub SSGeom
{
    my($pResults) = @_; ;
    my $result="";
    
    if ($$pResults{'SSGeom-BOOL'} ne "")
    {
        if    ($$pResults{'SSGeom-BOOL'} eq 'OK')  {$result='0';}
        elsif ($$pResults{'SSGeom-BOOL'} eq 'BAD') {$result='1';}    
    }
    else
    {  
        $result = '?';
    }
    
    return($result);
}
#-------------------------------------------------------------------------
sub Voids
{
    my($pResults) = @_;
    my $result1 = ""; 
    my $result2 = ""; my $result2Ref; my $Ntotal=0;
    my $result3 = ""; my $result3Ref; my $Mtotal=0;
    
    if (($$pResults{'Voids-NATIVE-LARGEST'} ne "") && ($$pResults{'Voids-MUTANT-LARGEST'} ne ""))
    {
        $result1 = $$pResults{'Voids-NATIVE-LARGEST'} - $$pResults{'Voids-MUTANT-LARGEST'}; 
        
    }
    else
    {
        $result1 = '0';  # shall i change this to zero? or set it to ?
    } 
    
    if ($$pResults{'Voids-NATIVE'} ne "")
    {
        $result2Ref = $$pResults{'Voids-NATIVE'}; # this is an array 
        $Ntotal = scalar (@$result2Ref);
        
        $result2 = readArrayRef($result2Ref,$Ntotal);

    }
    else 
    {
        $result2 = '0,0,0,0,0,0,0,0,0,0';
    } 
    
    if ($$pResults{'Voids-MUTANT'} ne "")
    {
        $result3Ref = $$pResults{'Voids-MUTANT'}; # this is an array
        $Mtotal = scalar (@$result3Ref);
        
        $result3 = readArrayRef($result3Ref,$Mtotal);
        
        }
    else
    {
        $result3 =  '0,0,0,0,0,0,0,0,0,0';
    }
    
    return($result1, $result2, $result3);        
}
#-------------------------------------------------------------------------
sub readArrayRef
{
my($result2Ref, $Ntotal) = @_;
my $result;
if    ($Ntotal == 10){$result ="$result2Ref->[0],$result2Ref->[1],$result2Ref->[2],$result2Ref->[3],$result2Ref->[4],$result2Ref->[5],$result2Ref->[6],$result2Ref->[7],$result2Ref->[8],$result2Ref->[9]";}
        elsif ($Ntotal == 9) {$result = "$result2Ref->[0],$result2Ref->[1],$result2Ref->[2],$result2Ref->[3],$result2Ref->[4],$result2Ref->[5],$result2Ref->[6],$result2Ref->[7],$result2Ref->[8],0";}
        elsif ($Ntotal == 8) {$result = "$result2Ref->[0],$result2Ref->[1],$result2Ref->[2],$result2Ref->[3],$result2Ref->[4],$result2Ref->[5],$result2Ref->[6],$result2Ref->[7],0,0"; }
        elsif ($Ntotal == 7) {$result ="$result2Ref->[0],$result2Ref->[1],$result2Ref->[2],$result2Ref->[3],$result2Ref->[4],$result2Ref->[5],$result2Ref->[6],0,0,0"; }
        elsif ($Ntotal == 6) {$result = "$result2Ref->[0],$result2Ref->[1],$result2Ref->[2],$result2Ref->[3],$result2Ref->[4],$result2Ref->[5],0,0,0,0";}
        elsif ($Ntotal == 5) {$result = "$result2Ref->[0],$result2Ref->[1],$result2Ref->[2],$result2Ref->[3],$result2Ref->[4],0,0,0,0,0";}
        elsif ($Ntotal == 4) {$result = "$result2Ref->[0],$result2Ref->[1],$result2Ref->[2],$result2Ref->[3],0,0,0,0,0,0";}
        elsif ($Ntotal == 3) {$result = "$result2Ref->[0],$result2Ref->[1],$result2Ref->[2],0,0,0,0,0,0,0"; }
        elsif ($Ntotal == 2) {$result = "$result2Ref->[0],$result2Ref->[1],0,0,0,0,0,0,0,0";  }
        else                 {$result = "$result2Ref->[0],0,0,0,0,0,0,0,0,0"; }
      return ($result);
}       
#-------------------------------------------------------------------------
sub Clash
{
    my($pResults) = @_;
    my $result="";
    
    if ($$pResults{'Clash-ENERGY'} ne "") 
    {
        $result = $$pResults{'Clash-ENERGY'}
    }
    else
    {
        $result = '?';
    }
    return($result);
}
#-------------------------------------------------------------------------
sub Glycine
{
    my($pResults) = @_;
    my $result="";
    
    if (($$pResults{'Glycine-NATIVE'} eq "GLY") && ($$pResults{'Glycine-MUTANT'} ne "GLY"))
    {
        if(($$pResults{'Glycine-NATIVE-ENERGY'} ne "") && ($$pResults{'Glycine-MUTANT-ENERGY'} ne ""))
        {
            $result = $$pResults{'Glycine-NATIVE-ENERGY'} - $$pResults{'Glycine-MUTANT-ENERGY'}; 
        }
        else 
        {
            $result = '?';
        }
    }
    else
    {
        #if not mutant to Gly we set the diffrent in energy to -10 or -100
        $result = '-100';
    }
    return($result);        
}
#-------------------------------------------------------------------------
sub Proline
{
    my($pResults) = @_;
    my $result = "";
    
    if (($$pResults{'Proline-MUTANT'} eq "PRO") && ($$pResults{'Proline-NATIVE'} ne "PRO"))
    {
        if (($$pResults{'Proline-MUTANT-ENERGY'} ne "") && ($$pResults{'Proline-NATIVE-ENERGY'} ne ""))
        {
            $result = $$pResults{'Proline-MUTANT-ENERGY'} - $$pResults{'Proline-NATIVE-ENERGY'}; 
        }
        else 
        {
            $result = '?';
        }
    }
    else
    {
        #if not mutant from Pro we set the diffrent in energy to -10 or -100
        $result = '-100';
    }
    return($result);        
}
#-------------------------------------------------------------------------
sub CisPro
{
    my($pResults) = @_;
    my $result = "";
    
    if ($$pResults{'CisPro-BOOL'} ne "")
    {
        if ($$pResults{'CisPro-BOOL'} eq 'OK')  {$result='0';}
        elsif ($$pResults{'CisPro-BOOL'} eq 'BAD') {$result='1';}    
    }
    else {$result = '?';}
    return($result);
}
#-------------------------------------------------------------------------