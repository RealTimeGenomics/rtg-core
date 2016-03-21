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
