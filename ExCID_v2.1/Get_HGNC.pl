#!/usr/bin/perl -w
use strict;
use LWP::Simple;
my $url = 'http://www.genenames.org/cgi-bin/download?'.
          'col=gd_app_sym&'.
          'col=gd_app_name&'.
          'col=gd_status&'.
          'col=gd_prev_sym&'.
          'col=gd_aliases&'.
          'col=gd_name_aliases&'.
          'col=gd_pub_chrom_map&'.
          'col=gd_pub_acc_ids&'.
          'col=gd_pub_ensembl_id&'.
          'col=gd_pub_refseq_ids&'.
          'col=gd_ccds_ids&'.
          'col=gd_vega_ids&'.
          'col=md_mim_id&'.
          'col=md_ucsc_id&'.
          'status=Approved&'.
          'status_opt=2&'.
          'where=&'.
          'order_by=gd_app_sym_sort&'.
          'format=text&'.
          'limit=&'.
          'submit=submit';
getprint($url);
