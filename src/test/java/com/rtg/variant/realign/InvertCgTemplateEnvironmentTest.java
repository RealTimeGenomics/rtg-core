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
package com.rtg.variant.realign;

import static com.rtg.util.StringUtils.LS;

import com.rtg.mode.DNA;
import com.rtg.util.TestUtils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.machine.MachineType;
import com.rtg.variant.VariantAlignmentRecord;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantParamsBuilder;
import com.rtg.variant.bayes.complex.ComplexTemplate;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import junit.framework.TestCase;

/**
 */
public class InvertCgTemplateEnvironmentTest extends TestCase {

  static final String TEMPLATE = ""
      //1234567890
      + "AAAAAAAAAA"
      + "CCCCCCCCCC"
      + "GGGGGGGGGG"
      + "TTTTTTTTTT"
      + "AAAAACCCCC"
      + "GGGGGT"
      ;
  static final String READ = ""
      + "ATTTTTTTTTTGGGGGGGGGGCCCCCCCCCCAAAA"
      ;
  static final String QUALITY = "\"%%%%%%%%%%``````````%%%%%%%%%%$`&$";

  public void test1() {
    final byte[] template = DNA.stringDNAtoByte(TEMPLATE);
    final byte[] read = DNA.stringDNAtoByte(READ);
    final SAMRecord sam = new SAMRecord(new SAMFileHeader());
    sam.setReadBases(READ.getBytes());
    sam.setAlignmentStart(4);
    sam.setCigarString("35=");
    sam.setFirstOfPairFlag(true);
    sam.setReadPairedFlag(true);
    sam.setReadNegativeStrandFlag(false);
    sam.setBaseQualityString(QUALITY);
    final VariantParams params = new VariantParamsBuilder().create();
    assertEquals(35, read.length);
    final AlignmentEnvironment temEnv = new AlignmentEnvironmentGenomeSubstitution(sam.getAlignmentStart() - 1, template.length, new ComplexTemplate(template, "", 20, 20), new byte[] {});
    final Environment env = new InvertCgTemplateEnvironment(
      new EnvironmentCombined(new AlignmentEnvironmentRead(new VariantAlignmentRecord(sam), params, MachineType.ILLUMINA_PE), sam.getAlignmentStart() - 1, 5, temEnv),
      MachineType.COMPLETE_GENOMICS);

    Exam.integrity(env);
    //    final String exp = ""
    //      + InvertCgTemplateEnvironment.class.getSimpleName() + " mReadStart=3 mMaxShift=5 mStart=20 mEnd=20 mReplacePosition=20 mDelta=0" + LS
    //      + "mSamEnv:" + AlignmentEnvironment.class.getSimpleName() + " read=[1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1] quality=[0.7943, 0.3981, 0.3981, 0.3981, 0.3981, 0.3981, 0.3981, 0.3981, 0.3981, 0.3981, 0.3981, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.3981, 0.3981, 0.3981, 0.3981, 0.3981, 0.3981, 0.3981, 0.3981, 0.3981, 0.3981, 0.5012, 0.0000, 0.3162, 0.5012]" + LS
    //      + "mTemplate:[56]" + LS
    //      + "[0] 1, 1, 1, 1, 1, 1, 1, 1, 1, 1" + LS
    //      + "[10] 2, 2, 2, 2, 2, 2, 2, 2, 2, 2" + LS
    //      + "[20] 3, 3, 3, 3, 3, 3, 3, 3, 3, 3" + LS
    //      + "[30] 4, 4, 4, 4, 4, 4, 4, 4, 4, 4" + LS
    //      + "[40] 1, 1, 1, 1, 1, 2, 2, 2, 2, 2" + LS
    //      + "[50] 3, 3, 3, 3, 3, 4" + LS
    //      + "mCot:" + ComplexTemplate.class.getSimpleName() + " mStart=20 mEnd=20" + LS
    //      + "mTemplate:[56]" + LS
    //      + "[0] 1, 1, 1, 1, 1, 1, 1, 1, 1, 1" + LS
    //      + "[10] 2, 2, 2, 2, 2, 2, 2, 2, 2, 2" + LS
    //      + "[20] 3, 3, 3, 3, 3, 3, 3, 3, 3, 3" + LS
    //      + "[30] 4, 4, 4, 4, 4, 4, 4, 4, 4, 4" + LS
    //      + "[40] 1, 1, 1, 1, 1, 2, 2, 2, 2, 2" + LS
    //      + "[50] 3, 3, 3, 3, 3, 4" + LS
    //      + "mReplaceBytes:[0]" + LS
    //      + "mReplaceString:" + LS
    //      + "" + LS
    //      + "mReplace:[0]" + LS
    //      + "template after replace:" + LS
    //      + "1, 1, 1, 1, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, " + LS
    //      ;
    final String x1 = ""
        + "template after replace:" + LS
        + "AAAAATTTTTTTTTTGGGGGGGGGGCCCCCCCCCCAAAAAAAAAANNNNNNNNNN"
        ;
    final String actual = env.toString();
    //System.err.println(actual);
    TestUtils.containsAll(actual, x1);
    assertEquals(35, env.readLength());
    assertEquals(41, env.absoluteTemplatePosition(0));
    assertEquals(37, env.absoluteTemplatePosition(4));
    assertEquals(7, env.absoluteTemplatePosition(34));


    assertEquals(1, env.template(-1));
    assertEquals(1, env.template(0));
    assertEquals(1, env.template(1));
    assertEquals(4, env.template(2));
    assertEquals(4, env.template(7));
    assertEquals(1, env.template(36));
    assertEquals(1, env.template(37));
  }

}
