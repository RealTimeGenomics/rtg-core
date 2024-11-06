/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
    return diversity.append(table).toString();
  }

  @Override
  public void generateReport() {
  }
}
