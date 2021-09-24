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

package com.rtg.simulation;

import java.io.File;
import java.io.IOException;

import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 * Test Class
 */
public class MutatedReferenceReadNameParserTest extends TestCase {

  @Override
  public void setUp() {
    Diagnostic.setLogStream();
  }

  public void testIndelCompensation() throws IOException {
    final File temp = FileHelper.createTempDirectory();
    try {
      final File mutations = FileHelper.resourceToFile("com/rtg/simulation/resources/mutations.vcf", new File(temp, "mutations.vcf"));
      checkArthur(mutations);
      checkMartha(mutations);
      checkPaul(mutations);
      checkSusan(mutations);
    } finally {
      assertTrue(FileHelper.deleteAll(temp));
    }
  }

  private void checkArthur(File mutations) throws IOException {
    final SimulatedReadNameParser arthurParser = new MutatedReferenceReadNameParser(new NewReadNameParser(), MutatedSampleOffsets.getOffsets(mutations, "Arthur"));

    check(arthurParser, "read111 0/3/simulatedSequence1_0/1506/R/101.", "simulatedSequence1", 1506);
    check(arthurParser, "read112 0/3/simulatedSequence1_0/1507/R/101.", "simulatedSequence1", 1507); //Closest position rounding up
    check(arthurParser, "read113 0/3/simulatedSequence1_0/1508/R/101.", "simulatedSequence1", 1507);

    check(arthurParser, "read211 0/3/simulatedSequence1_0/9778/R/101.", "simulatedSequence1", 9777);
    check(arthurParser, "read212 0/3/simulatedSequence1_0/9779/R/101.", "simulatedSequence1", 9780);

    check(arthurParser, "read120 0/3/simulatedSequence1_1/40/R/101.", "simulatedSequence1", 40);
    check(arthurParser, "read120 0/3/simulatedSequence1_1/3502/R/101.", "simulatedSequence1", 3502);
    check(arthurParser, "read121 0/3/simulatedSequence1_1/3503/R/101.", "simulatedSequence1", 3503);
    check(arthurParser, "read122 0/3/simulatedSequence1_1/3504/R/101.", "simulatedSequence1", 3503); //Closest position rounding up
    check(arthurParser, "read123 0/3/simulatedSequence1_1/3505/R/101.", "simulatedSequence1", 3504); //Closest position rounding up
    check(arthurParser, "read124 0/3/simulatedSequence1_1/3506/R/101.", "simulatedSequence1", 3504); //Closest position rounding up
    check(arthurParser, "read125 0/3/simulatedSequence1_1/3507/R/101.", "simulatedSequence1", 3504);
    check(arthurParser, "read126 0/3/simulatedSequence1_1/3508/R/101.", "simulatedSequence1", 3505);
    check(arthurParser, "read126 0/3/simulatedSequence1_1/5000/R/101.", "simulatedSequence1", 4997);

    check(arthurParser, "read1506098 0/3/simulatedSequence2_1/62/F/101.", "simulatedSequence2", 62);

    check(arthurParser, "read311 0/3/simulatedSequence13/182/R/101.", "simulatedSequence13", 182);
    check(arthurParser, "read312 0/3/simulatedSequence13/183/R/101.", "simulatedSequence13", 183); //Closest position rounding up
    check(arthurParser, "read313 0/3/simulatedSequence13/184/R/101.", "simulatedSequence13", 183);

    check(arthurParser, "read411 0/3/simulatedSequence14/393/R/101.", "simulatedSequence14", 393);
    check(arthurParser, "read412 0/3/simulatedSequence14/394/R/101.", "simulatedSequence14", 393); //Closest position rounding up
    check(arthurParser, "read413 0/3/simulatedSequence14/395/R/101.", "simulatedSequence14", 394); //Closest position rounding up
    check(arthurParser, "read414 0/3/simulatedSequence14/396/R/101.", "simulatedSequence14", 394);

    check(arthurParser, "read511 0/3/simulatedSequence14/1151/R/101.", "simulatedSequence14", 1149);
    check(arthurParser, "read512 0/3/simulatedSequence14/1152/R/101.", "simulatedSequence14", 1152);

    check(arthurParser, "read2197622 0/26/simulatedSequence15/54/F/101.", "simulatedSequence15", 54);
    check(arthurParser, "read611 0/26/simulatedSequence15/528/F/101.", "simulatedSequence15", 528);
    check(arthurParser, "read612 0/26/simulatedSequence15/529/F/101.", "simulatedSequence15", 528); //Closest position rounding up
    check(arthurParser, "read613 0/26/simulatedSequence15/530/F/101.", "simulatedSequence15", 528); //Closest position rounding up
    check(arthurParser, "read614 0/26/simulatedSequence15/531/F/101.", "simulatedSequence15", 528); //Closest position rounding up
    check(arthurParser, "read615 0/26/simulatedSequence15/532/F/101.", "simulatedSequence15", 529); //Closest position rounding up
    check(arthurParser, "read616 0/26/simulatedSequence15/533/F/101.", "simulatedSequence15", 529); //Closest position rounding up
    check(arthurParser, "read617 0/26/simulatedSequence15/534/F/101.", "simulatedSequence15", 529); //Closest position rounding up
    check(arthurParser, "read618 0/26/simulatedSequence15/535/F/101.", "simulatedSequence15", 529);
  }

  private void checkMartha(File mutations) throws IOException {
    final SimulatedReadNameParser marthaParser = new MutatedReferenceReadNameParser(new NewReadNameParser(), MutatedSampleOffsets.getOffsets(mutations, "Martha"));

    check(marthaParser, "read211 0/3/simulatedSequence1_0/9800/R/101.", "simulatedSequence1", 9800);

    check(marthaParser, "read120 0/3/simulatedSequence1_1/1506/F/101.", "simulatedSequence1", 1506);
    check(marthaParser, "read121 0/3/simulatedSequence1_1/1507/F/101.", "simulatedSequence1", 1507); //Closest position rounding up
    check(marthaParser, "read122 0/3/simulatedSequence1_1/1508/F/101.", "simulatedSequence1", 1507);

    check(marthaParser, "read123 0/3/simulatedSequence1_1/3533/F/101.", "simulatedSequence1", 3532);
    check(marthaParser, "read124 0/3/simulatedSequence1_1/3534/F/101.", "simulatedSequence1", 3537);

    check(marthaParser, "read1506098 0/3/simulatedSequence2_1/62/F/101.", "simulatedSequence2", 62);

    check(marthaParser, "read311 0/3/simulatedSequence13_0/1258/R/101.", "simulatedSequence13", 1258);
    check(marthaParser, "read312 0/3/simulatedSequence13_0/1259/R/101.", "simulatedSequence13", 1260);

    check(marthaParser, "read313 0/3/simulatedSequence13_1/9000/R/101.", "simulatedSequence13", 9000);

    check(marthaParser, "read2197622 0/26/simulatedSequence15/54/F/101.", "simulatedSequence15", 54);
    check(marthaParser, "read411 0/26/simulatedSequence15/532/F/101.", "simulatedSequence15", 532);
    check(marthaParser, "read412 0/26/simulatedSequence15/533/F/101.", "simulatedSequence15", 534);
  }

  private void checkPaul(File mutations) throws IOException {
    final SimulatedReadNameParser paulParser = new MutatedReferenceReadNameParser(new NewReadNameParser(), MutatedSampleOffsets.getOffsets(mutations, "Paul"));

    check(paulParser, "read111 0/3/simulatedSequence1_0/1506/R/101.", "simulatedSequence1", 1506);
    check(paulParser, "read112 0/3/simulatedSequence1_0/1507/R/101.", "simulatedSequence1", 1507); //Closest position rounding up
    check(paulParser, "read113 0/3/simulatedSequence1_0/1508/R/101.", "simulatedSequence1", 1507);

    check(paulParser, "read121 0/3/simulatedSequence1_1/1506/F/101.", "simulatedSequence1", 1506);
    check(paulParser, "read122 0/3/simulatedSequence1_1/1507/F/101.", "simulatedSequence1", 1507); //Closest position rounding up
    check(paulParser, "read123 0/3/simulatedSequence1_1/1508/F/101.", "simulatedSequence1", 1507);

    check(paulParser, "read211 0/3/simulatedSequence1_0/3535/R/101.", "simulatedSequence1", 3534);

    check(paulParser, "read221 0/3/simulatedSequence1_1/3533/R/101.", "simulatedSequence1", 3532);
    check(paulParser, "read222 0/3/simulatedSequence1_1/3534/R/101.", "simulatedSequence1", 3537);

    check(paulParser, "read311 0/3/simulatedSequence1_0/9778/R/101.", "simulatedSequence1", 9777);
    check(paulParser, "read311 0/3/simulatedSequence1_0/9779/R/101.", "simulatedSequence1", 9780);

    check(paulParser, "read321 0/3/simulatedSequence1_1/9774/R/101.", "simulatedSequence1", 9777);
    check(paulParser, "read322 0/3/simulatedSequence1_1/9775/R/101.", "simulatedSequence1", 9778);
    check(paulParser, "read322 0/3/simulatedSequence1_1/9776/R/101.", "simulatedSequence1", 9779);

    check(paulParser, "read1506098 0/3/simulatedSequence2_1/62/F/101.", "simulatedSequence2", 62);

    check(paulParser, "read411 0/3/simulatedSequence13/1258/R/101.", "simulatedSequence13", 1258);
    check(paulParser, "read412 0/3/simulatedSequence13/1259/R/101.", "simulatedSequence13", 1260);

    check(paulParser, "read511 0/3/simulatedSequence14/393/R/101.", "simulatedSequence14", 393);
    check(paulParser, "read512 0/3/simulatedSequence14/394/R/101.", "simulatedSequence14", 393); //Closest position rounding up
    check(paulParser, "read513 0/3/simulatedSequence14/395/R/101.", "simulatedSequence14", 394); //Closest position rounding up
    check(paulParser, "read514 0/3/simulatedSequence14/396/R/101.", "simulatedSequence14", 394);

    check(paulParser, "read521 0/3/simulatedSequence14/1151/R/101.", "simulatedSequence14", 1149);
    check(paulParser, "read522 0/3/simulatedSequence14/1152/R/101.", "simulatedSequence14", 1152);

    check(paulParser, "read2197622 0/26/simulatedSequence15/54/F/101.", "simulatedSequence15", 54);
    check(paulParser, "read2197622 0/26/simulatedSequence15/54/F/101.", "simulatedSequence15", 54);
    check(paulParser, "read611 0/26/simulatedSequence15/528/F/101.", "simulatedSequence15", 528);
    check(paulParser, "read612 0/26/simulatedSequence15/529/F/101.", "simulatedSequence15", 528); //Closest position rounding up
    check(paulParser, "read613 0/26/simulatedSequence15/530/F/101.", "simulatedSequence15", 528); //Closest position rounding up
    check(paulParser, "read614 0/26/simulatedSequence15/531/F/101.", "simulatedSequence15", 528); //Closest position rounding up
    check(paulParser, "read615 0/26/simulatedSequence15/532/F/101.", "simulatedSequence15", 529); //Closest position rounding up
    check(paulParser, "read616 0/26/simulatedSequence15/533/F/101.", "simulatedSequence15", 529); //Closest position rounding up
    check(paulParser, "read617 0/26/simulatedSequence15/534/F/101.", "simulatedSequence15", 529); //Closest position rounding up
    check(paulParser, "read618 0/26/simulatedSequence15/535/F/101.", "simulatedSequence15", 529);
  }

  private void checkSusan(File mutations) throws IOException {
    final SimulatedReadNameParser susanParser = new MutatedReferenceReadNameParser(new NewReadNameParser(), MutatedSampleOffsets.getOffsets(mutations, "Susan"));

    check(susanParser, "read111 0/3/simulatedSequence1_0/1510/R/101.", "simulatedSequence1", 1510);

    check(susanParser, "read120 0/3/simulatedSequence1_1/1506/F/101.", "simulatedSequence1", 1506);
    check(susanParser, "read121 0/3/simulatedSequence1_1/1507/F/101.", "simulatedSequence1", 1507); //Closest position rounding up
    check(susanParser, "read122 0/3/simulatedSequence1_1/1508/F/101.", "simulatedSequence1", 1507);

    check(susanParser, "read211 0/3/simulatedSequence1_0/3503/R/101.", "simulatedSequence1", 3503);
    check(susanParser, "read212 0/3/simulatedSequence1_0/3504/R/101.", "simulatedSequence1", 3503); //Closest position rounding up
    check(susanParser, "read213 0/3/simulatedSequence1_0/3505/R/101.", "simulatedSequence1", 3504); //Closest position rounding up
    check(susanParser, "read214 0/3/simulatedSequence1_0/3506/R/101.", "simulatedSequence1", 3504); //Closest position rounding up
    check(susanParser, "read215 0/3/simulatedSequence1_0/3507/R/101.", "simulatedSequence1", 3504);

    check(susanParser, "read311 0/3/simulatedSequence1_1/3533/F/101.", "simulatedSequence1", 3532);
    check(susanParser, "read312 0/3/simulatedSequence1_1/3534/F/101.", "simulatedSequence1", 3537);

    check(susanParser, "read411 0/3/simulatedSequence1_0/9000/F/101.", "simulatedSequence1", 8997);
    check(susanParser, "read412 0/3/simulatedSequence1_1/9000/F/101.", "simulatedSequence1", 9003);

    check(susanParser, "read1506098 0/3/simulatedSequence2_1/62/F/101.", "simulatedSequence2", 62);

    check(susanParser, "read511 0/3/simulatedSequence13_0/9000/R/101.", "simulatedSequence13", 9000);

    check(susanParser, "read521 0/3/simulatedSequence13_1/182/R/101.", "simulatedSequence13", 182);
    check(susanParser, "read522 0/3/simulatedSequence13_1/183/R/101.", "simulatedSequence13", 183); //Closest position rounding up
    check(susanParser, "read523 0/3/simulatedSequence13_1/184/R/101.", "simulatedSequence13", 183);

    check(susanParser, "read531 0/3/simulatedSequence13_1/1259/R/101.", "simulatedSequence13", 1258);
    check(susanParser, "read532 0/3/simulatedSequence13_1/1260/R/101.", "simulatedSequence13", 1260);

    check(susanParser, "read541 0/3/simulatedSequence13_1/9000/R/101.", "simulatedSequence13", 9000);

    check(susanParser, "read2197622 0/26/simulatedSequence15/54/F/101.", "simulatedSequence15", 54);
    check(susanParser, "read611 0/26/simulatedSequence15/532/F/101.", "simulatedSequence15", 532);
    check(susanParser, "read612 0/26/simulatedSequence15/533/F/101.", "simulatedSequence15", 534);
  }

  private void check(SimulatedReadNameParser parser, String read, String expectedTemplateName, long expectedPosition) {
    parser.setReadInfo(read, -1);
    assertEquals(expectedTemplateName, parser.templateName());
    assertEquals(expectedPosition, parser.templatePosition());
  }
}
