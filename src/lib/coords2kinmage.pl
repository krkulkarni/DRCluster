#!/usr/bin/perl
# coords2kinmage.pl
use strict;

open F, $ARGV[0];

my %P;

while ( <F> )
{
   chomp;
   my @a = split /\t/;
   my ($id) = $a[0] =~ /(PF\d+)/i;
   push @{ $P{$id} }, sprintf("%.3f", $a[1]) . ' ' . sprintf("%.3f", $a[2]) . ' ' . sprintf("%.3f", $a[3]);
}
close F;

my $ids = keys %P;
$ids += 1 if $ids % 2 != 0;

my @hsv_gradient;
for (my $i=0; $i<360; $i+=360/($ids/2))
{
   unshift @hsv_gradient, "$i 100 100";
   push @hsv_gradient, "$i 100 50";
}

print "\@kinemage 1\n";
print "\@whitebackground\n";
 
my $i = 0;
foreach my $id (keys %P)
{
   print "\@hsvcolor {$id} $hsv_gradient[$i]\n";
   $i++;
}

foreach my $id (keys %P)
{
   print "\@balllist {$id} color= $id radius= 0.25 nohighlight\n"; 
   
   foreach my $point ( @{ $P{$id} } )
   {
      print "$point\n";
   }
}
