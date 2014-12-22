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
package com.rtg.variant.bayes.multisample;

import static com.rtg.sam.SharedSamConstants.SAM_CANCER;
import static com.rtg.util.StringUtils.LS;

import java.io.File;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.sam.SamFilterOptions;
import com.rtg.tabix.TabixIndexer;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

/**
 */
public class MultisampleCliTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new MultisampleCli();
  }

  public void testInitParams() {
    checkHelp("multisnp [OPTION]... -o DIR -r FILE -t SDF FILE+",
              "SAM/BAM format files containing mapped reads",
              "directory for output",
              "input-list-file",
              "average coverage determined from calibration files",
              "of the reference genome the reads"
    );
  }

  private static final String CALIBRATION = ""
    + "@mnp:RG1\t0\t5\t1\t1\t0\t0\t4\t1" + LS
    + "@nh:RG1\t0\t6" + LS
    + "@covar\treadgroup\tbasequality\tsequence\tequal\tdiff\tins\tdel" + LS
    + "RG1\t63\tg1\t7\t41\t0\t0" + LS;

  public void testValidator() throws Exception {
    try (final TestDirectory tmpDir = new TestDirectory()) {
      final File tmpFile = FileUtils.stringToFile("original-derived TEST cancer contamination=0.13", new File(tmpDir, "relations.relations"));
      String err = checkMainInitBadFlags("-o", tmpDir.getPath(), "-t", tmpDir.getPath(), "-r", tmpFile.getPath());
      assertTrue(err, err.contains("The directory \"" + tmpDir.getPath() + "\" already exists"));

      final File outDir = new File(tmpDir, "out");
      err = checkMainInitBadFlags("-o", outDir.getPath(), "-t", tmpDir.getPath(), "-r", tmpFile.getPath());
      assertTrue(err, err.contains("No input files specified"));

      final File template = new File(tmpDir, "template");
      ReaderTestUtils.getDNADir(">g1\nacgtacgtacgtacgtacgt", template);

      err = checkMainInitBadFlags("-o", outDir.getPath(), "-t", template.getPath(), "-r", tmpFile.getPath());
      assertTrue(err, err.contains("No input files specified"));

      final File in = new File(tmpDir, "alignments.sam.gz");
      FileHelper.stringToGzFile(SAM_CANCER, in);
      new TabixIndexer(in, new File(tmpDir, "alignments.sam.gz.tbi")).saveSamIndex();

      err = checkMainInitBadFlags("-o", outDir.getPath(), "-t", template.getPath(), "-r", tmpFile.getPath(), "-T", "-1", in.getPath());
      assertTrue(err, err.contains("flag \"--threads\" has invalid value \"-1\""));

      err = checkMainInitBadFlags("-o", outDir.getPath(), "-t", template.getPath(), "-r", tmpFile.getPath(), "--" + SamFilterOptions.MAX_HITS_FLAG, "-1", in.getPath());
      assertTrue(err, err.contains("flag \"--max-hits\" has invalid value \"-1\""));

      final File cal = new File(tmpDir, "alignments.sam.gz.calibration");
      FileUtils.stringToFile(CALIBRATION, cal);

      err = checkMainInitBadFlags("-o", outDir.getPath(), "-t", template.getPath(), "-r", tmpFile.getPath(), in.getPath(), in.getPath(), "--filter-depth", "100", "--filter-depth-multiplier", "3.0");
      assertTrue(err, err.contains("Only one of --filter-depth or --filter-depth-multiplier can be set"));

      err = checkMainInitBadFlags("-o", outDir.getPath(), "-t", template.getPath(), "-r", tmpFile.getPath(), in.getPath(), in.getPath(), "--max-coverage", "100", "--max-coverage-multiplier", "3.0");
      assertTrue(err, err.contains("Only one of --max-coverage or --max-coverage-multiplier can be set"));

      err = checkMainInitBadFlags("-o", outDir.getPath(), "-t", template.getPath(), "-r", tmpFile.getPath(), in.getPath(), in.getPath(), "--filter-depth", "100", "--max-coverage-multiplier", "3.0");
      assertTrue(err, err.contains("Cannot mix ratio based and fixed coverage thresholding"));

      err = checkMainInitBadFlags("-o", outDir.getPath(), "-t", template.getPath(), "-r", tmpFile.getPath(), in.getPath(), in.getPath(), "--max-coverage", "100", "--filter-depth-multiplier", "3.0");
      assertTrue(err, err.contains("Cannot mix ratio based and fixed coverage thresholding"));

      err = checkMainInitBadFlags("-o", outDir.getPath(), "-t", template.getPath(), "-r", tmpFile.getPath(), in.getPath(), in.getPath(), "--Xlookahead", "-1");
      assertTrue(err, err.contains("Invalid lookahead value, minimum is 2"));

      err = checkMainInitBadFlags("-o", outDir.getPath(), "-t", template.getPath(), "-r", tmpFile.getPath(), in.getPath(), in.getPath(), "--Xchunking", "-1");
      assertTrue(err, err.contains("Invalid chunk size, minimum chunking size is 1000"));

      final File bed = new File(tmpDir, "some.bed");
      err = checkMainInitBadFlags("-o", outDir.getPath(), "-t", template.getPath(), "-r", tmpFile.getPath(), in.getPath(), in.getPath(), "--filter-bed", bed.getPath());
      assertTrue(err, err.contains("File not found: \"The specified file, \"" + bed.getPath() + "\", does not exist"));
      assertTrue(bed.createNewFile());
      checkMainInitOk("-o", outDir.getPath(), "-t", template.getPath(), "-r", tmpFile.getPath(), in.getPath(), in.getPath(), "--filter-bed", bed.getPath());
      FileHelper.deleteAll(outDir);

      err = checkMainInitWarn("-o", outDir.getPath(), "-t", template.getPath(), "-r", tmpFile.getPath(), in.getPath(), in.getPath(), "--filter-depth-multiplier", "10.0", "--max-coverage-multiplier", "5.0");
      assertTrue(err, err.contains("--filter-depth-multiplier is set higher than --max-coverage-multiplier, consider increasing --max-coverage-multiplier"));

      FileHelper.deleteAll(outDir);


      err = checkMainInitWarn("-o", outDir.getPath(), "-t", template.getPath(), "-r", tmpFile.getPath(), in.getPath(), in.getPath(), "--filter-depth", "100", "--max-coverage", "10");
      assertTrue(err, err.contains("--filter-depth is set higher than --max-coverage, consider increasing --max-coverage"));


      FileHelper.deleteAll(outDir);

      checkMainInitOk("-o", outDir.getPath(), "-t", template.getPath(), "-r", tmpFile.getPath(), in.getPath(), in.getPath());
    }
  }
}
