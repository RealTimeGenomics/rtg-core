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
package com.rtg.simulation;

import com.rtg.util.StringUtils;



/**
 * Parses names generated by the <code>DWGsim</code> program.
 * <code>http://sourceforge.net/apps/mediawiki/dnaa/index.php?title=Whole_Genome_Simulation</code>
 *
 * <code>@&lt;#1&gt;_&lt;#2&gt;_&lt;#3&gt;_&lt;#4&gt;_&lt;#5&gt;_&lt;#6&gt;_&lt;#7&gt;_&lt;#8&gt;:&lt;#9&gt;:&lt;#10&gt;_&lt;#11&gt;:&lt;#12&gt;:&lt;#13&gt;_&lt;#14&gt;</code>
 * #1 contig name (chromosome name)
 * #2 start end 1 (zero-based)
 * #3 start end 2 (zero-based)
 * #4 strand end 1 (0 - forward, 1 - reverse)
 * #5 strand end 2 (0 - forward, 1 - reverse)
 * #6 random read end 1 (0 - from the mutated reference, 1 - random)
 * #7 random read end 2 (0 - from the mutated reference, 1 - random)
 * #8 number of sequencing errors end 1 (colour errors for colorspace)
 * #9 number of SNPs end 1
 * #10 number of indels end 1
 * #11 number of sequencing errors end 2 (colour errors for colorspace)
 * #12 number of SNPs end 2
 * #13 number of indels end 2
 * #14 read number (unique within a given contig/chromosome)
 *
 * e.g. <code>@CFTR.3.70s_282_148_0_0_0_0_4:0:0_4:0:0_10c/1</code>
 *
 */
public class DwgsimReadNameParser implements SimulatedReadNameParser {

  static final char READ_SIM_SEPARATOR = '_';
  static final int EXPECT_FIELDS = 10;

  /** Index of the template sequence id field. */
  public static final int TEMPLATE_SEQ_NAME = 0;
  /** Index of the template position field 1. */
  public static final int TEMPLATE_POS_ID_1 = 1;
  /** Index of the template position field 2. */
  public static final int TEMPLATE_POS_ID_2 = 2;
  /** Index of strand direction field 1 */
  public static final int STRAND_DIR_1 = 3;
  /** Index of strand direction field 2 */
  public static final int STRAND_DIR_2 = 4;

  /** Index of number of SNPs field 1 */
  public static final int NUM_ERRORS_1 = 7;
  /** Index of number of indels field 1 */
  public static final int NUM_ERRORS_2 = 8;

  /** Index of the read label field. */
  public static final int READ_LABEL = 9;

  /**
   * @param readName an example read name
   * @return true if this read name looks like something this parser deals with
   */
  public static boolean looksOk(String readName) {
    if (readName.indexOf(READ_SIM_SEPARATOR) != -1) {
      final String[] f = StringUtils.split(readName, READ_SIM_SEPARATOR);
      if (f.length == EXPECT_FIELDS) {
        return true;
      }
    }
    return false;
  }



  private String[] mFields = null;
//  private int mSubs;
//  private int mIns;
//  private int mDel;
  private int mNumMismatches;
  private boolean mFirstArm;
  private long mReadNum;

  @Override
  public boolean setReadInfo(String readName, int readLength) {
    mFields = StringUtils.split(readName, READ_SIM_SEPARATOR);
    if (mFields.length != EXPECT_FIELDS) {
      return false;
    }

      mFirstArm = mFields[READ_LABEL].charAt(mFields[READ_LABEL].length() - 1) == '1';
      final int slashpos = mFields[READ_LABEL].indexOf('/');
      if (slashpos == -1) {
        return false;
      }
      mReadNum = Long.parseLong(mFields[READ_LABEL].substring(0, slashpos), 16);

      final String[] errors = StringUtils.split(mFields[mFirstArm ? NUM_ERRORS_1 : NUM_ERRORS_2], ':');
      if (errors.length != 3) {
        return false;
      }
      mNumMismatches = Integer.parseInt(errors[0]);
    return true;
  }

  @Override
  public boolean isChimera() {
    return false;
  }

  @Override
  public boolean isDuplicate() {
    return false;
  }

  @Override
  public long templatePosition() {
    if (mFirstArm) {
      return Long.parseLong(mFields[TEMPLATE_POS_ID_1]);
    } else {
      return Long.parseLong(mFields[TEMPLATE_POS_ID_2]);
    }
  }

  @Override
  public String templateName() {
    return mFields[TEMPLATE_SEQ_NAME];
  }

  @Override
  public int templateSet() {
    return 0; // This assumes the generator is operating in a haploid mode.
  }

  @Override
  public boolean forwardFrame() {
    if (mFirstArm) {
      return Integer.parseInt(mFields[STRAND_DIR_1]) == 0;
    } else {
      return Integer.parseInt(mFields[STRAND_DIR_2]) == 0;
    }
  }

  @Override
  public int readId() {
    return (int) mReadNum; // Ick
  }

  @Override
  public String readName() {
    return String.valueOf(mReadNum);
  }

  @Override
  public int substitutions() {
    return 0;
  }

  @Override
  public int insertions() {
    return mNumMismatches;  //this is stupid, but for 'worst score' purposes, probably a good idea.
  }

  @Override
  public int deletions() {
    return 0; //this is rolled into insertions. Should work out fine.
  }

  @Override
  public int numMismatches() {
    return mNumMismatches;
  }
}
