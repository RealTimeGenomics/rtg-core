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
package com.rtg.sam;

import java.io.File;
import java.io.IOException;

import com.rtg.util.InvalidParamsException;
import com.rtg.util.TestUtils;
import com.rtg.util.cli.CFlags;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 */
public class SamCommandHelperTest extends TestCase {

  private CFlags mFlags = null;

  @Override
  public void setUp() {
    mFlags = new CFlags();
    SamCommandHelper.initSamRg(mFlags);
  }

  @Override
  public void tearDown() {
    mFlags = null;
  }


  public void testRG() throws IOException {
    final File tmpDir = FileUtils.createTempDir("tmp", "sdjf");
    try {
      mFlags.setFlags("--sam-rg", tmpDir.getAbsolutePath());
      assertFalse(SamCommandHelper.validateSamRg(mFlags));
      assertTrue(mFlags.getParseMessage().contains("for --sam-rg is a directory, must be a file"));

      final File f = new File(tmpDir, "samrgfile");
      FileUtils.stringToFile("", f);

      mFlags.setFlags("--sam-rg", f.getAbsolutePath());
      assertTrue(SamCommandHelper.validateSamRg(mFlags));
    } finally {
      FileHelper.deleteAll(tmpDir);
    }
  }

    public void testSamRGErrors() throws IOException, InvalidParamsException {
   final File outer = FileUtils.createTempDir("rammap", "end2end");
   try {

       final File header = new File(outer, "header");
       FileUtils.stringToFile("", header);

       final MemoryPrintStream mps = new MemoryPrintStream();
       Diagnostic.setLogStream(mps.printStream());
       try {
         try {
           SamCommandHelper.validateAndCreateSamRG(header.toString(), SamCommandHelper.ReadGroupStrictness.REQUIRED);
           fail();
         } catch (final InvalidParamsException ipe) {
           assertTrue(ipe.getMessage().contains("file \"" + header.getPath()));
//           assertTrue(mps.toString().contains("No read group information present in the input file \"" + header.getPath() + "\", please provide file with single read group line"));
         }
         assertNull(SamCommandHelper.validateAndCreateSamRG(header.toString(), SamCommandHelper.ReadGroupStrictness.OPTIONAL));

         final File header2 = new File(outer, "header2");
         FileUtils.stringToFile("@RG\tID:L23\tSM:NA123" + "\n" + "@RG\tID:L43\tSM:NA123", header2);

         final MemoryPrintStream mps2 = new MemoryPrintStream();
         Diagnostic.setLogStream(mps2.printStream());
         try {
           SamCommandHelper.validateAndCreateSamRG(header2.toString(), SamCommandHelper.ReadGroupStrictness.REQUIRED);
           fail();
         } catch (final InvalidParamsException ipe) {
           assertTrue(ipe.getMessage().contains("file \"" + header2.getPath()));
//           assertTrue(mps2.toString().contains("Multiple read group information present in the input file \"" + header2.getPath() + "\", please provide file with single read group line"));
         }
         try {
           SamCommandHelper.validateAndCreateSamRG(header2.toString(), SamCommandHelper.ReadGroupStrictness.AT_MOST_ONE);
           fail();
         } catch (final InvalidParamsException ipe) {
           assertTrue(ipe.getMessage().contains("file \"" + header2.getPath()));
//           assertTrue(mps2.toString().contains("Multiple read group information present in the input file \"" + header2.getPath() + "\", please provide file with single read group line"));
         }
         assertNull(SamCommandHelper.validateAndCreateSamRG(header2.toString(), SamCommandHelper.ReadGroupStrictness.OPTIONAL));
       } finally {
         Diagnostic.setLogStream();
       }
   } finally {
     assertTrue(FileHelper.deleteAll(outer));
   }
  }

  public void testStrictnessEnum() {
    TestUtils.testEnum(SamCommandHelper.ReadGroupStrictness.class, "[REQUIRED, OPTIONAL, AT_MOST_ONE]");
  }

}
