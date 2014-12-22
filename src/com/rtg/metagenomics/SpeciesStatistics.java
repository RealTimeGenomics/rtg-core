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

package com.rtg.metagenomics;

import java.io.File;

import com.rtg.launcher.AbstractStatistics;
import com.rtg.util.StringUtils;
import com.rtg.util.TextTable;

/**
 */
public class SpeciesStatistics extends AbstractStatistics {

  private static final String FORMAT = "%1.5g";
  double mRichness = 0;
  double mShannon = 0;
  double mPielou = 0;
  double mInvSimpson = 0;

  /**
   * @param outputDirectory The base output directory to generate statistics and reports in. May be null if no statistics or reports are to be generated.
   */
  public SpeciesStatistics(File outputDirectory) {
    super(outputDirectory);
  }

  @Override
  protected String getStatistics() {
    final StringBuilder diversity = new StringBuilder("Diversity metrics").append(StringUtils.LS);
    final TextTable table = new TextTable();
    table.addRow("Richness:", String.format(FORMAT, mRichness));
    table.addRow("Shannon:", String.format(FORMAT, mShannon));
    table.addRow("Pielou:", String.format(FORMAT, mPielou));
    table.addRow("Inverse Simpson:", String.format(FORMAT, mInvSimpson));
    return diversity.append(table.toString()).toString();
  }

  @Override
  public void generateReport() {
  }
}
