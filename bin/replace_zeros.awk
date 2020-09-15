#!/bin/awk -f
<<<<<<< HEAD
=======

>>>>>>> 3ed2bddd5686fcd937992d9f55274d29ea7fd703
BEGIN{FS=OFS="\t"}

{
  delete a;
  nzero=0;

  for (i=4;i<=NF;i++)
  {
    if ($i==0) a[++nzero]=i;
    if ($i!=0 && $i!=1 && $i!="NA")
    {
      print;
      next;
    }
  }
  for (i=1;i<=nzero;i++) $a[i]=0.001;
  print;
}
