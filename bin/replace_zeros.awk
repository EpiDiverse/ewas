#!/bin/awk -f

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
