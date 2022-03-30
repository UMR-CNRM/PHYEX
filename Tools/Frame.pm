package Tools::Frame;

use strict;
use base qw (Exporter);
our @EXPORT = qw (frame);

sub frame
{
  my ($t, $width, $border) = @_;
  my $len = length ($t);

  $border = 1 unless (defined ($border));

  my $S = ' ';

  my ($C, $V, $H) = ('*', '|', '-');

  ($C, $V, $H) = (' ', ' ', ' ')
    unless ($border);

  my $df = 3;
  $width ||= $len + 2 * $df;

  my $line1 = $C . ($H x ($width-2)) . $C;
  my $line2 = $V . ($S x ($width-2)) . $V;


  my $TEXT = '';

  $TEXT .= "$line1\n";
  for (1 .. ($df-1)/2)
    {
      $TEXT .= "$line2\n";
    }

  die ("Cannot frame text: `$t'\n")
    if ($width - 2 * $df <= 0);
  
  while ($t)
    {
      my $s = substr ($t, 0, $width - 2 * $df, '');

      my $i = 0;
      while (length ($s) < $width - 2 * $df)
        {
          if ($i % 2)
            {
              $s = " $s";
            }
          else
            {
              $s = "$s ";
            }
          $i++;
        }
      my $linet = $V . ($S x ($df-1)) . $s .  ($S x ($df-1)) . $V;
      $TEXT .= "$linet\n";
    }

  for (1 .. ($df-1)/2)
    {
      $TEXT .= "$line2\n";
    }
  $TEXT .= "$line1\n";
}


1;
