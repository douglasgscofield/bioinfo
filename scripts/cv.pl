# This program accepts the output of 'show-coords -lcdTH' and displays the mapping query(s) to the reference(s)

##########################################
#	Programmers: Alan Twaddle, University of Pittsburgh
#				Center for Vaccine Research
#					Ghedin Lab
#	Date: February 1, 2012


#!/usr/bin/perl -w

use strict;
use warnings;
use Switch;
use Getopt::Long;
use Bio::Graphics::Panel;
use Bio::Graphics::Feature;

my $cc = 0;
my $count = 0;
my $ref_len = 0;
my $ftr = 'Bio::Graphics::Feature';
my $temp_color = "";
my @seq_array = ();
my $temp_name = "";
my @read_array = ();
my $read = "";
my $seg = "";
my $whole_seq;
my $ref = "";
my $query = "";
my $temp_query = "";
#my %options=();
#getopts("ql",\%options);


#################### Command Line Arguments ####################

my $HELP_INFO = q~
  USAGE: perl cv.pl [options] <coords-output>

  DESCRIPTION:
    Displays mummer/nucmer mappings and prints PNG to stdout via the "show-coords -lcdTH" 
    output. The output file is <coords-output>.png. Command line options allow for
    different features when producing the image. Within the code there are comment
    blocks which give examples of how to further customize it i.e. grouping based
    on name similarity, color, etc. This program is aimed at providing full
    customization within the limits of the Bio-Graphics module.
    Current version supported is Bio-Graphics-2.25 

  OPTIONS:
    -q		    Groups segments of the same query sequence together into the same
		    track, i.e. transcript (default off)
    -l		    Adds query sequence label (default off)
    -g
    --grid	    Toggles vertical grid bars in background (default off)
    --grid_color    Set color of grid bars (default light cyan)
    --grid_major    Set color of major grid bars (default cyan)
    --glyph	    Selects type of feature to draw (default generic).
    -c
    --connector	    Chooses the type of connector when features are grouped (default dashed)
    -d		    Toggles display of feature orientation (default off)
    -w
    --width	    Sets width of image in pixels (default 600)
    --pad_left	    Pad left side of image (white space) in pixels (default 10)
    --pad_right	    Pad right side of image (white space) in pixels (default 10)
    -f		    Flip  each feature so that largest coordinate is leftmost coordinate
    -s
    --spacing	    Spacing between each track (in pixels) (default 2)
    --colored	    Invoke coloring scheme (default Reference=light_blue, Queries=red)
		    When invoked the default is each feature is assigned a different
		    color. When used with '-q' the entire feature is assigned a uniq color
    --key	    Draw a key at the bottom of panel, between each track, left/right of 
		    each track, or none (default none)
    --font	    Changes the font of the query labels. See GD::Text for font options.
    --scale_font    Changes font for scale label.
    --key_font	    Changes key font labels.
    --key_color	    Changes background color for key.
    -h
    --help	    Prints this help message

~;

my $help;
my $q_group;
my $label;
my $grid = 0;
my $grid_color = "cyan";
my $grid_major = "light cyan";
my $glyph = "generic";
my $dir;
my $connector = "dashed";
my $width = 600;
my $pl = 10;
my $pr = 10;
my $flip = 0;
my $space = 2;
my $col;
my $colored;
my $key = "none";
my $font = "gdSmallFont";
my $key_color = "white";
my $scale_font = "gdSmallFont";
my $key_font = "gdSmallFont";

GetOptions( 
	'q' => \$q_group,
	'l' => \$label,
	'g|grid' => \$grid,
	'grid_color:s' => \$grid_color,
	'glyph:s' => \$glyph,
	'd' => \$dir,
	'c|connector:s' => \$connector,
	'width|w:i' => \$width,
	'pad_left:i' => \$pl,
	'pad_righti:i' => \$pr,
	'f' => \$flip,
	'spacing|s:i' => \$space,
	'collapse' =>\$col,
	'colored' => \$colored,
	'key:s' => \$key,
	'grid_major:s' => \$grid_major,
	'font:s' => \$font,
	'key_color:s' => \$key_color,
	'scale_font:s' => \$scale_font,
	'key_font:s' => \$key_font,
	'h|help' => \$help
	
);

if(defined $help){
	print $HELP_INFO;
	exit;
}

#################### Command Line Arguments ####################



################### Color Generator ###################
my @colors;

if(defined $colored){

	for (my $i = 0; $i < 255; $i++) {

		my ($rand,$x);
		my @hex;

		for ($x = 0; $x < 3; $x++) {
		$rand = rand(255);
		$hex[$x] = sprintf ("%x", $rand);

			if ($rand < 9) {

				$hex[$x] = "0" . $hex[$x];

			}

			if ($rand > 9 && $rand < 16) {

				$hex[$x] = "0" . $hex[$x];
			
			}
		}

		$colors[$i] = "\#" . $hex[0] . $hex[1] . $hex[2];
	}
}


################### Color Generator ###################

################### Create Feature ###################

sub newFeature{

	chomp;
	my @fa = split(/\t/,$_);

	my $nf = $ftr->new(-start=>$fa[0],
			   -end=>$fa[1],
			   -strand=>$fa[12],
			   -name=> $fa[14]);

	return $nf;

}

################### Create Feature ###################

open IN, "<$ARGV[0]";
open OUT,">$ARGV[0].png";
my @in = <IN>;
close IN;

# Obtains width of largest reference so that it can be displayed properly

foreach(@in){

	my @line = split("\t",$_);
	
	if($count == 0){

		$ref_len = $line[7];
		$count = 1;

	}elsif($line[7] > $ref_len){

		$ref_len = $line[7];

	}

}


# Creates Panel
my $panel = Bio::Graphics::Panel->new(
	-length=>$ref_len+$pl, 
	-width=>$width, 
	-spacing=>$space,
	-grid=>$grid,
	-gridcolor=>$grid_color, 
	-key_style=>"bottom",
	-key_color=>$key_color, 
	-pad_left=>$pl, 
	-pad_rigt=>$pr,
	-flip=> $flip,
	-key_style=>"$key",
	-gridmajorcolor=>"$grid_major"
);


# Retrieves color based on the name of feature

sub getColor{

	$cc += 1;

	if(defined $colored){
=begin
		if($temp_color eq ""){
		
			$temp_color = $_[0];
			return $colors[$cc];

		}
		elsif($temp_color eq $_[0]){

			return $colors[$cc];

		}
		else{

			$temp_color = $_[0];
			$cc += 1;
			return $colors[$cc];

		}

=cut

# Use the switch case block for grouping features together based on similar names

		switch($_[0]){
			case /014/	{ return "red"; next;}
			case /002/	{ return "blue"; next;}
			case /006/	{ return "green"; next;}
			case /007/	{ return "yellow"; next;}
			case /023/	{ return "orange"; next;}
		#	case /Imm/	{ return $colors[5]; next;}
			else		{ return "purple";	}
		}

	}else{
		return "red";
	}

}

#my $panel = Bio::Graphics::Panel->new(-length=>2900, -width=>2900, -spacing=>2, -grid=>1, -gridcolor=>'lightgrey', -key_style=>"bottom", -pad_left=>50, -pad_rigt=>50);


my $arrow = $ftr->new(-start=>0, -end=>$ref_len);

$panel->add_track($arrow, -tick=>2, -glyph=>'arrow', -font=>$scale_font);

foreach(@in){
	
	@seq_array = split("\t",$_);

	$ref = $seq_array[13];
	
	if($temp_name eq ""){

                $temp_name = $seq_array[13];
		$whole_seq = $ftr->new(-start=>"0", -end=>$seq_array[7]);
                $panel->add_track($whole_seq, -label=>$ref, -key=>"$ref Reference", -font=>$font);

        }
	
	if($ref ne $temp_name){

		$read = $ftr->new(-segments=>\@read_array, -name=>$temp_query);

		$panel->add_track($read,-bgcolor=> getColor($temp_query),
					-label=>
						sub{if(defined $label){
							return 1;
						}},
					 -connector=>$connector,
					 -font => $font
					# -key=>$temp_query
				);


		@read_array = ();

		$whole_seq = $ftr->new(-start=>"0", -end=>$seq_array[7]);
		$panel->add_track($whole_seq, -label=>$ref, -key => "$ref Reference", -font=>$font);
		$temp_name = $ref;

	}

#	$seg = $ftr->new(-start=>$seq_array[0], -end=>$seq_array[1], -source=>$ref, -strand=>$seq_array[12]);
#
	$seg = newFeature($_);

	if(defined $q_group){
		if($temp_query eq ""){

			$temp_query = $seq_array[14];
		
		}

		if($temp_query ne $seq_array[14]){

			$read = $ftr->new(-segments=>\@read_array, -name=>$temp_query);

	                $panel->add_track(
				$read,
				-bgcolor=> getColor($temp_query),
				-connector=>$connector, 
				-label=>
					sub{if(defined $label){
						return 1;
					}},
				 -glyph=>$glyph,
				 -font => $font
				#-key=>$temp_query
			);

			@read_array = ();


		}

		$temp_query = $seq_array[14];

		push(@read_array, $seg);

	}else{

		$panel->add_track($seg,-bgcolor=> getColor($seq_array[14]),
					-glyph=>$glyph, 
					-label=>
						sub{if(defined $label){
							return 1; 
						}},
					-font => $font
				#	-key => $temp_query
				);
	}

}

$read = $ftr->new(-segments=>\@read_array, -name=>$temp_query);

$panel->add_track($read,-bgcolor=>getColor($temp_query),
			-connector=>$connector, 
			-label=>
				sub{if(defined $label){
					return 1;
				}},
			-glyph=>$glyph,
			-font => $font
		#	-key => $temp_query
		);




# Use this to create key features if you are using grouped features of some sort

#$panel->add_track(generic => [],
#		-key=>"INS002",
#		-bgcolor=>"blue");
$panel->add_track(generic => [],
                -key=>"INS006",
                -bgcolor=>"red",
		-font=>$key_font);
$panel->add_track(generic => [],
                -key=>"INS007",
                -bgcolor=>"green",
		-font=>$key_font);
$panel->add_track(generic => [],
                -key=>"INS014",
                -bgcolor=>"yellow",
		-font=>$key_font);
$panel->add_track(generic => [],
                -key=>"INS023",
                -bgcolor=>"orange",
		-font=>$key_font);



print OUT $panel->png;
close OUT;
