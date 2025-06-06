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

package com.rtg.variant.sv;

import static com.rtg.util.StringUtils.LS;

import java.io.File;
import java.io.IOException;

import com.rtg.launcher.SequenceParams;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.TestUtils;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.Flag;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;
import com.rtg.variant.sv.SvCliUtils.SvValidator;
import com.rtg.variant.sv.SvParamsTest.MockSvParams;
import com.rtg.variant.sv.SvParamsTest.MockSvParams.MockSvParamsBuilder;

import junit.framework.TestCase;

/**
 */
public class SvCliUtilsTest extends TestCase {

  public void testInitFlags() {
    final CFlags flags = new CFlags();
    SvCliUtils.initCommonFlags(flags);
    SvCliUtils.initRelabelFlag(flags);
    TestUtils.containsAll(flags.getCompactFlagUsage()
        , "[OPTION]... -o DIR -t SDF -r FILE FILE+"
        , "[OPTION]... -o DIR -t SDF -I FILE -R FILE"
        );
    TestUtils.containsAllUnwrapped(flags.getUsageString()
        , "-I,", "--input-list-file=FILE", "file containing a list of SAM/BAM format files (1 per line) containing mapped reads"
        , "-o,", "--output=DIR", "directory for output"
        , "--readgroup-labels=FILE", "file containing read group relabel mappings (1 per line)"
        , "-r,", "--readgroup-stats=FILE", "text file containing read group stats. May be specified 0 or more times"
        , "-R,", "--readgroup-stats-list-file=FILE", "file containing list of read group stats files (1 per line)"
        , "-t,", "--template=SDF", "reference genome"
        , "FILE+", "SAM/BAM format files containing mapped reads. May be specified 0 or more times"
        , "-m,", "--max-as-mated=INT", "if set, ignore mated SAM records with an alignment score (AS attribute) that exceeds this value"
        , "-u,", "--max-as-unmated=INT", "if set, ignore unmated SAM records with an alignment score (AS attribute) that exceeds this value"
        , "--region=REGION", "if set, only process SAM records within the specified range. "
        , "-h,", "--help", "print help on command-line flag usage"
        , "-Z,", "--no-gzip", "do not gzip the output"
        , "-T,", "--threads=INT", "number of threads (Default is the number of available cores)"
        );
    final Flag<?> regionFlag = flags.getFlag("region");
    assertNull(regionFlag.getChar());
  }

  private static class MockSvValidator extends SvValidator { }

  public void testPopulateCommonParamsAndValidator() throws IOException, InvalidParamsException {
    Diagnostic.setLogStream();
    final CFlags flags = new CFlags();
    flags.setName("blah");
    SvCliUtils.initCommonFlags(flags);
    flags.setValidator(new MockSvValidator());
    flags.setInvalidFlagHandler(null);
    final File tempDir = FileUtils.createTempDir("discordantTool", "makeParamsTest");
    try {
      final File f = FileUtils.stringToFile("RG1\t200000\t56893\t1991255\t69693925\t55590\t22261070\t9022416526\t55590\t22261070\t9022416526\t5\t55590\t236\t1067", new File(tempDir, "blah.txt"));
      final File gen = ReaderTestUtils.getDNADir(">g1" + LS + "aaatcgactggtcagctagg" + LS, tempDir);
      final File out = new File(tempDir, "blah");
      assertTrue(flags.setFlags("--template", gen.getPath(), "--output", out.getPath(), "-r", f.getPath(), f.getPath(), "--threads", "4", "--no-gzip"));

      MockSvParamsBuilder builder = new MockSvParamsBuilder();
      SvCliUtils.populateCommonParams(builder, SequenceParams.builder(), flags);
      MockSvParams params = builder.create();
      assertEquals("blah", params.name());
      assertEquals(3, params.ioThreads());
      assertEquals(out.getPath(), params.directory().getPath());
      assertFalse(params.outputParams().isCompressed());

      assertTrue(flags.setFlags("--template", gen.getPath(), "--output", out.getPath(), "-r", f.getPath(), f.getPath(), "--threads", "4"));
      builder = new MockSvParamsBuilder();
      SvCliUtils.populateCommonParams(builder, SequenceParams.builder(), flags);
      params = builder.create();
      assertTrue(params.outputParams().isCompressed());
      final File listFile = FileUtils.stringToFile(f.getPath(), new File(tempDir, "list.txt"));
      assertTrue(flags.setFlags("--template", gen.getPath(), "--output", out.getPath(), "-R", listFile.getPath(), "-I", listFile.getPath()));
      assertFalse(flags.setFlags("--template", gen.getPath(), "--output", f.getPath(), "-r", f.getPath(), f.getPath()));
      assertFalse(flags.setFlags("--template", f.getPath(), "--output", out.getPath(), "-r", f.getPath(), f.getPath()));
      assertFalse(flags.setFlags("--template", gen.getPath(), "--output", out.getPath(), "-R", out.getPath(), f.getPath()));
      assertFalse(flags.setFlags("--template", gen.getPath(), "--output", out.getPath(), "-r", f.getPath(), "-I", out.getPath()));
      assertFalse(flags.setFlags("--template", gen.getPath(), "--output", out.getPath(), "-r", f.getPath(), f.getPath(), "--threads", "0"));
      assertFalse(flags.setFlags("--template", gen.getPath(), "--output", out.getPath(), "-r", f.getPath(), f.getPath(), "--region", "blah:12-11"));
    } finally {
      assertTrue(FileHelper.deleteAll(tempDir));
    }
  }

  public void testInvalidStats() throws Exception {

    final CFlags flags = new CFlags();
    SvCliUtils.initCommonFlags(flags);
    final File tempDir = FileUtils.createTempDir("discordantTool", "makeParamsTest");
    try {
      final File f = FileUtils.stringToFile("RG1\t200000\t56893\t1991255\t69693925\t55590\t22261070\t9022416526\t55590\t22261070\t9022416526\t5\t55590\t236\t1067", new File(tempDir, "blah.txt"));
      flags.setFlags("--" + SvCliUtils.RG_STATS_FILE, f.getPath(), "-t", "irrelevant", "-o", "alsoirrelevant");

      final MockSvParamsBuilder builder = new MockSvParamsBuilder();
      SvCliUtils.populateReadGroupStats(builder, flags);

      try {
        FileUtils.stringToFile("RG1\t0\t56893\t1991255\t69693925\t55590\t22261070\t9022416526\t55590\t22261070\t9022416526\t5\t55590\t236\t1067", new File(tempDir, "blah.txt"));
        SvCliUtils.populateReadGroupStats(builder, flags);
        fail();
      } catch (final NoTalkbackSlimException ntse) {
      }
      try {
        FileUtils.stringToFile("RG1\t200000\t0\t1991255\t69693925\t55590\t22261070\t9022416526\t55590\t22261070\t9022416526\t5\t55590\t236\t1067", new File(tempDir, "blah.txt"));
        SvCliUtils.populateReadGroupStats(builder, flags);
        fail();
      } catch (final NoTalkbackSlimException ntse) {
      }
      try {
        FileUtils.stringToFile("RG1\t200000\t56893\t1991255\t69693925\t0\t22261070\t9022416526\t55590\t22261070\t9022416526\t5\t55590\t236\t1067", new File(tempDir, "blah.txt"));
        SvCliUtils.populateReadGroupStats(builder, flags);
        fail();
      } catch (final NoTalkbackSlimException ntse) {
      }
      try {
        FileUtils.stringToFile("RG1\t200000\t56893\t1991255\t69693925\t55590\t22261070\t9022416526\t0\t22261070\t9022416526\t5\t55590\t236\t1067", new File(tempDir, "blah.txt"));
        SvCliUtils.populateReadGroupStats(builder, flags);
        fail();
      } catch (final NoTalkbackSlimException ntse) {
      }
      try {
        FileUtils.stringToFile("RG1\t200000\t56893\t1991255\t69693925\t55590\t22261070\t9022416526\t55590\t22261070\t9022416526\t5\t-2\t236\t1067", new File(tempDir, "blah.txt"));
        SvCliUtils.populateReadGroupStats(builder, flags);
        fail();
      } catch (final NoTalkbackSlimException ntse) {
      }
      try {
        FileUtils.stringToFile("RG1\t200000\t56893\t1991255\t69693925\t55590\t22261070\t9022416526\t55590\t22261070\t9022416526\t5\t55590\t-2\t1067", new File(tempDir, "blah.txt"));
        SvCliUtils.populateReadGroupStats(builder, flags);
        fail();
      } catch (final NoTalkbackSlimException ntse) {
      }
      try {
        FileUtils.stringToFile("RG1\t200000\t56893\t1991255\t69693925\t55590\t22261070\t9022416526\t55590\t22261070\t9022416526\t5\t55590\t236\t-2", new File(tempDir, "blah.txt"));
        SvCliUtils.populateReadGroupStats(builder, flags);
        fail();
      } catch (final NoTalkbackSlimException ntse) {
      }
    } finally {
      assertTrue(FileHelper.deleteAll(tempDir));
    }
  }
}
