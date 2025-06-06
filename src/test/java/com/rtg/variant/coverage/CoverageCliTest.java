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
package com.rtg.variant.coverage;

import static com.rtg.sam.SharedSamConstants.REF_SEQS_M;
import static com.rtg.sam.SharedSamConstants.SAM_M;
import static com.rtg.util.StringUtils.LS;

import java.io.File;
import java.io.IOException;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.tabix.TabixIndexer;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.TestUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

/**
 */
public class CoverageCliTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new CoverageCli();
  }

  private void checkParamsError(final String[] args, final String... exp) {
    TestUtils.containsAllUnwrapped(checkHandleFlagsErr(args), exp);
  }

  public void testErrorFlags() {
    checkParamsError(new String[] {}, "Error: You must provide a value for -o DIR");
  }

  public void testInvalidParams() throws IOException, InvalidParamsException {
    final File temp = FileUtils.createTempDir("test", "temp");
    try {
      final File out = new File(temp, "out");
      final File template = ReaderTestUtils.getDNADir(REF_SEQS_M, temp);
      final File sam = new File(temp, "alignments.sam.gz");
      FileHelper.stringToGzFile(SAM_M, sam);
      checkParamsError(new String[] {"-o", sam.getPath(), "-t", template.getPath()}, "The directory \"" + sam.getPath() + "\" already exists. Please remove it first or choose a different directory.");
      checkParamsError(new String[] {"-o", out.getPath(), "-t", out.getPath()}, "The specified SDF, \"" + out.getPath() + "\", does not exist.");
      checkParamsError(new String[] {"-o", out.getPath(), "-t", sam.getPath()}, "The specified file, \"" + sam.getPath() + "\", is not an SDF.");
      checkParamsError(new String[] {"-o", out.getPath(), "-t", template.getPath(), "-s", "-1"}, "--smoothing must be in the range [0, 8192]");
      checkParamsError(new String[] {"-o", out.getPath(), "-t", template.getPath(), "--input-list-file", out.getPath()});
      checkParamsError(new String[] {"-o", out.getPath(), "-t", template.getPath(), sam.getPath(), "-T", "0"});
      checkParamsError(new String[] {"-o", out.getPath(), "-t", template.getPath(), sam.getPath(), "-m", "-1"});
    } finally {
      assertTrue(FileHelper.deleteAll(temp));
    }
  }

  public void testValidParams() throws Exception {
    final File temp = FileUtils.createTempDir("test", "temp");
    try {
      final File out = new File(temp, "out");
      final File template = ReaderTestUtils.getDNADir(REF_SEQS_M, temp);
      final File sam = new File(temp, "alignments.sam.gz");
      FileHelper.stringToGzFile(SAM_M, sam);
      final File inputList = new File(temp, "list.txt");
      FileUtils.stringToFile(sam.getPath() + LS, inputList);
      new TabixIndexer(sam, new File(temp, "alignments.sam.gz.tbi")).saveSamIndex();
      checkMainInitOk("-o", out.getPath(), "-t", template.getPath(), "-I", inputList.getPath(), "--region", "g1");
      assertTrue(out.exists());
      assertTrue(new File(out, "coverage.bed.gz").exists());
      final String outputBed = FileHelper.gzFileToString(new File(out, "coverage.bed.gz"));
      assertFalse(outputBed, outputBed.contains("g2"));
    } finally {
      assertTrue(FileHelper.deleteAll(temp));
    }
  }

  public void testValidParamsNoRef() throws Exception {
    final File temp = FileUtils.createTempDir("test", "temp");
    try {
      final File out = new File(temp, "out");
      final File sam = new File(temp, "alignments.sam.gz");
      FileHelper.stringToGzFile(SAM_M, sam);
      final File inputList = new File(temp, "list.txt");
      FileUtils.stringToFile(sam.getPath() + LS, inputList);
      new TabixIndexer(sam, new File(temp, "alignments.sam.gz.tbi")).saveSamIndex();
      checkMainInitWarn("-o", out.getPath(), "-I", inputList.getPath(), "--region", "g1");
      assertTrue(out.exists());
      assertTrue(new File(out, "coverage.bed.gz").exists());
      final String outputBed = FileHelper.gzFileToString(new File(out, "coverage.bed.gz"));
      assertFalse(outputBed, outputBed.contains("g2"));
    } finally {
      assertTrue(FileHelper.deleteAll(temp));
    }
  }

  public void testInitParams() {
    checkHelp("coverage [OPTION]... -o DIR FILE+"
            , "Measures and reports coverage depth of read alignments across a reference."
            , "SAM/BAM format files containing mapped reads"
            , "-o,", "--output=DIR", "directory for output"
            , "-t,", "--template=SDF", "reference genome"
            , "-I,", "--input-list-file=FILE", "file containing a list of SAM/BAM format files (1 per line) containing mapped reads"
            , "\\0, --exclude-mated", "exclude all mated SAM records"
            , "\\0, --exclude-unmated", "exclude all unmated SAM records"
            , "-m,", "--max-as-mated=INT", "if set, ignore mated SAM records with an alignment score (AS attribute) that exceeds this value"
            , "-u,", "--max-as-unmated=INT", "if set, ignore unmated SAM records with an alignment score (AS attribute) that exceeds this value"
            , "-c,", "--max-hits=INT", "if set, ignore SAM records with an alignment count that exceeds this value"
            , "\\0, --region=REGION", "if set, only process SAM records within the specified range. "
            , "-s,", "--smoothing", "smooth with this number of neighboring values (0 means no smoothing) (Default is 50)"
            , "-Z,", "--no-gzip", "do not gzip the output"
            , "--per-base", "output per-base counts in TSV format"
            , "--per-region", "output BED/BEDGRAPH entries per-region"
            , "-T,", "--threads=INT", "number of threads (Default is the number of available cores)"
            );
    checkExtendedHelp("coverage [OPTION]... -o DIR FILE+"
            , "--Xerror-rates", "report statistics about sequencer error rates"
            , "--Xcoverage-threshold=INT"
            , "--Xbinarize-bed"
            , "--Xcallability"
            );
  }
}
