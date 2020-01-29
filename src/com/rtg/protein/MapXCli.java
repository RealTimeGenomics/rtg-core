/*
 * Copyright (c) 2014. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */
package com.rtg.protein;

/**
 * Module for doing protein matching using translated reads
 */
public class MapXCli extends MapProteinCli {

  @Override
  public String moduleName() {
    return "mapx";
  }

  @Override
  public String description() {
    return "searches translated nucleotide read data sets against a protein database";
  }

  @Override
  protected boolean translated() {
    return true;
  }

  @Override
  protected String queryLabel() {
    return "read";
  }
}
