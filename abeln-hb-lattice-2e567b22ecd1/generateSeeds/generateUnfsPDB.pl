#!/usr/bin/perl

use strict;

use POSIX qw(ceil floor);

my @ABC=('A' .. 'Z');

### FILES ###
my $fn_faa = $ARGV[0];
my $num_seeds=$ARGV[1];
my $fn_AA="AtoAaa.txt";
my $output_dir="output/";

my $NUM_RES=0;

### GLOBALS ###

my %AtoAaa=();
my @SEQUENCE=();

### PROGRAM ###

readAA();
readFaa();
writePDB();
################


sub readAA{
    open(AA,$fn_AA) or die "can't open $fn_AA \n";
    while(my $line = <AA>){
	chomp $line;
	my ($A,$Aaa)=split(/\s/,$line);
	$AtoAaa{$A}=$Aaa;
    }
    close AA;
    print "read $fn_AA\n";
}

sub readFaa{
    #if($fn_faa !~ /\.faa/){ die "incorrect file type: $fn_faa\n"};
    open(FAA,$fn_faa) or die "can't open $fn_faa \n";
    while(my $line = <FAA>){
	chomp $line;
	if($line =~ /^>/){ next;}
	@SEQUENCE = split(//,$line);
    }
    $NUM_RES=$#SEQUENCE +1 ;
    print "$NUM_RES residues \n";
    close FAA;
    print "read $fn_faa\n";
}


sub writePDB{
    my @fnids = split('\/',$fn_faa);
    my $fnid = $fnids[$#fnids];
    print  $fnid . "\n";
    my ($name,$tmp)= split(/\./,$fnid);
    my $fn_out= $output_dir.'unfs'.$num_seeds.'_'.$name.'.pdb';
    open(OUT,">".$fn_out) or die "can't open $fn_out";
    my $n=0;
    #printf OUT "12345678901234567890123456789012345678901234567890123456789012345678\n";
    print OUT "MODEL\n";
    for(my $cn=0;$cn<$num_seeds;$cn++){ 
	$n=0;
	foreach my $A (@SEQUENCE){
	    my $Aaa =$AtoAaa{$A};
	    writeAA($Aaa,$n,$cn);
	    $n++;
	}
	print OUT "TER\n";
    }
    #print OUT "TER\n";
    close $fn_out;
    print "written to $fn_out\n";
}





sub writeAA{
    my ($AA,$resnum,$chainnum)=@_;
    if($resnum > 29){print "too many residues\n";}    

    my $startX=3;
    my $startY=3;
    my $startZ=3;

    my $atomNameMain = "CA";
    my $atomNameSide = "CB";
    my $Natom=$resnum;
    #my $atomName="O";
    my $aminoAcid=$AA;
    my $chain=$ABC[$chainnum];
    my $Nresidue=$resnum;
    
    my ($xcoord,$ycoord,$zcoord);
    my $sdir=1;
    if($chainnum < ($num_seeds/2)){
	$xcoord = $startX+3*($resnum);
	$ycoord = $startY ;
	$zcoord = $startZ +15*$chainnum;
	if($resnum % 2 ==0){$sdir=-1;}
    }else{
	$xcoord = $startX+3*($resnum);
	$ycoord = $startY +66;
	$zcoord = $startZ +15*($chainnum - floor(int($num_seeds)/int(2))); 
	if($resnum % 2 ==1){$sdir=-1;}
    }

    my $spinx =$xcoord;
    my $spiny = $ycoord +$sdir;
    my $spinz = $zcoord;

    my $inBeta =0;
    if($resnum !=0 && $resnum != ($NUM_RES -1)){
	$inBeta=1;
    }

###### BACKBONE
#1-6
    printf OUT "ATOM  " ;                         
#7-11
    printf OUT "% 5d", $Natom;                      
    printf OUT "";             
#13-16                   
    printf OUT "%4.4s",$atomNameMain;                    
    printf OUT "  ";  
#18-20  
    printf OUT "%3.3s",$aminoAcid;
    printf OUT " ";
#22
    printf OUT "%1.1s",$chain;
#23-26
    printf OUT "%4.4d", $Nresidue;
    printf OUT "    ";
#21-38
    printf OUT "%8.3f",$xcoord;
#39-46
    printf OUT "%8.3f",$ycoord;
#47-54
    printf OUT "%8.3f",$zcoord;	  
#55-60 occupancy
    printf OUT "%6.2f",1.0;
#61-66 b-factor
    if($inBeta){
	printf OUT "%6.2f",11.0;
    }else{
	printf OUT "%6.2f",22.0;
    }
    printf OUT "  \n";
###### SPIN / SIDE CHAIN
#1-6
    printf OUT "ATOM  " ;                         
#7-11
    printf OUT "% 5d", $Natom;                      
    printf OUT "";             
#13-16                   
    printf OUT "%4.4s",$atomNameSide;                    
    printf OUT "  ";  
#18-20  
    printf OUT "%3.3s",$aminoAcid;
    printf OUT " ";
#22
    printf OUT "%1.1s",$chain;
#23-26
    printf OUT "%4.4d", $Nresidue;
    printf OUT "    ";
#21-38
    printf OUT "%8.3f",$spinx;
#39-46
    printf OUT "%8.3f",$spiny;
#47-54
    printf OUT "%8.3f",$spinz;	  
#55-60 occupancy
    printf OUT "%6.2f",1.0;
#61-66 b-factor
    if($inBeta){
	printf OUT "%6.2f",11.0;
    }else{
	printf OUT "%6.2f",22.0;
    }
    printf OUT "  \n";
}






