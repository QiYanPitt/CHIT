#!/usr/bin/perl -w

$input1=$ARGV[0];
$output_file=$ARGV[1];

open INPUT1, "$input1";
open OUTPUT, ">>$output_file";
@raw_data1=<INPUT1>;

@broken0 = split(/\t/, $raw_data1[0]);
$start = $broken0[0];
$end = $broken0[1];
for ($i=0; $i<=$#raw_data1; $i++)
{  @broken1 = split(/\t/, $raw_data1[$i]);
   $j = $i+1;
   @broken2 = split(/\t/, $raw_data1[$j]);
   if($broken1[2] eq $broken2[2]){ 
     $start = $start . ';' . $broken2[0];
     $end = $end . ';'. $broken2[1];
   }
   if($broken1[2] ne $broken2[2]){
     print OUTPUT "$start $end $broken1[2]";
     $start = $broken2[0];
     $end = $broken2[1];
   }
}
