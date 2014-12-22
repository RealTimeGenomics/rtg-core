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

import static com.rtg.util.StringUtils.LS;
import static com.rtg.util.StringUtils.TAB;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;

import com.rtg.Slim;
import com.rtg.launcher.AbstractParamsCliTest;
import com.rtg.launcher.GlobalFlags;
import com.rtg.launcher.ParamsCli;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SdfId;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.sam.SharedSamConstants;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.test.FileHelper;

/**
 */
public class SpeciesCliTest extends AbstractParamsCliTest<SpeciesParams> {

  private static final String APP_NAME = "rtg species";

  @Override
  public void setUp() throws IOException {
    super.setUp();
    GlobalFlags.resetAccessedStatus();
    //comment this out if you want diagnostics from the test
    Diagnostic.setLogStream();
  }

  @Override
  public final void testApplicationName() {
    assertEquals(APP_NAME, new SpeciesCli().applicationName() + " " + new SpeciesCli().moduleName());
  }

  public void testHeaderExplicitly() {
    // When testing the header elsewhere, just use SpeciesTask.SPECIES_HEADER to minimize
    // code changes when column names change.
    assertEquals("abundance\tabundance-low\tabundance-high\tDNA-fraction\tDNA-fraction-low\tDNA-fraction-high\tconfidence\tcoverage-depth\tcoverage-breadth\treference-length\tmapped-reads\thas-reference\ttaxa-count\ttaxon-id\tparent-id\trank\ttaxonomy-name", SpeciesTask.SPECIES_HEADER);
  }

  public final void testInitFlags() {
    checkHelp(APP_NAME
        , "Calculates a species distribution from a metagenomic sample."
        , "-t,", "--genomes=SDF", "SDF containing the genomes"
        , "-I,", "--input-list-file=FILE", "file containing a list of SAM/BAM format files (1 per line) containing mapped reads"
        , "-o,", "--output=DIR", "directory for output"
        , "-r,", "--relabel-species-file=FILE", "file containing list of species name to reference name mappings (1 mapping per line format: [reference short name][tab][species])"
        , "FILE+", "SAM/BAM format files containing mapped reads. May be specified 0 or more times"
        , "--exclude-mated", "exclude all mated SAM records"
        , "--exclude-unmated", "exclude all unmated SAM records"
        , "--print-all", "print non present species in the output file"
        , "-c", "--min-confidence", "species below this confidence value will not be reported"
        , "-m,", "--max-as-mated=INT", "if set, ignore mated SAM records with an alignment score (AS attribute) that exceeds this value"
        , "-u,", "--max-as-unmated=INT", "if set, ignore unmated SAM records with an alignment score (AS attribute) that exceeds this value"
        , "-h,", "--help", "print help on command-line flag usage"
        );
    checkExtendedHelp(APP_NAME
        , "--Xhelp", "print help on extended command-line flag usage"
        , "--Xverbose", "turn on output of convergence information"
        , "--Xiterations=INT", "minimum number of iterations multiplied by the block size (Default is 1)"
        );
  }

  @Override
  protected ParamsCli<SpeciesParams> getParamsCli() {
    return new SpeciesCli();
  }

  public void testFileLoad() throws IOException, InvalidParamsException {
    final File dir = FileUtils.createTempDir("test", "speciesCli");
    try {
      final File t = new File(dir, "t");
      ReaderTestUtils.getReaderDNA(">a\nacgt", t, new SdfId(1L)).close();
      final File o = new File(dir, "o");
      final File list = new File(dir, "list");
      final File f1 = new File(dir, "f1");
      final File f2 = new File(dir, "f2");
      final File f3 = new File(dir, "f3");
      final File f4 = new File(dir, "f4");
      assertTrue(f1.createNewFile());
      assertTrue(f2.createNewFile());
      assertTrue(f3.createNewFile());
      assertTrue(f4.createNewFile());
      FileUtils.stringToFile(f1.getPath() + "\n" + f2.getPath() + "\n", list);
      SpeciesParams p = checkMakeParamsOut("-t", t.getPath(), "-o", o.getPath(), "-I", list.getPath());
      assertTrue(p.mapped().contains(f1));
      assertTrue(p.mapped().contains(f2));
      assertFalse(p.mapped().contains(f3));
      assertFalse(p.mapped().contains(f4));
      p = checkMakeParamsOut("-t", t.getPath(), "-o", o.getPath(), "-I", list.getPath(), f3.getPath(), f4.getPath());
      assertTrue(p.mapped().contains(f1));
      assertTrue(p.mapped().contains(f2));
      assertTrue(p.mapped().contains(f3));
      assertTrue(p.mapped().contains(f4));
      p = checkMakeParamsOut("-t", t.getPath(), "-o", o.getPath(), f3.getPath(), f4.getPath());
      assertFalse(p.mapped().contains(f1));
      assertFalse(p.mapped().contains(f2));
      assertTrue(p.mapped().contains(f3));
      assertTrue(p.mapped().contains(f4));

    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  private static final String SAM_HEADER = ""
      + "@HD" + TAB + "VN:1.0" + TAB + "SO:coordinate\n"
      + "@SQ" + TAB + "SN:a1" + TAB + "LN:4" + LS
      + "@SQ" + TAB + "SN:a2" + TAB + "LN:4" + LS
      + "@SQ" + TAB + "SN:b1" + TAB + "LN:4" + LS
      + "@SQ" + TAB + "SN:b2" + TAB + "LN:4" + LS
      + "@CO" + TAB + "TEMPLATE-SDF-ID:00000000-0000-0000-0000-000000000001";

  private static final String SAM = ""
      + SAM_HEADER + LS
      + "0" + TAB + "0" + TAB + "a1" + TAB + "1" + TAB + "255" + TAB + "4=" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ACGT" + TAB + "*" + LS
      + "1" + TAB + "0" + TAB + "a2" + TAB + "1" + TAB + "255" + TAB + "4=" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ACGT" + TAB + "*" + LS
      + "2" + TAB + "0" + TAB + "b1" + TAB + "1" + TAB + "255" + TAB + "4=" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ACGT" + TAB + "*" + LS
      + "3" + TAB + "0" + TAB + "b2" + TAB + "1" + TAB + "255" + TAB + "4=" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ACGT" + TAB + "*";

  static final String EXTRA_SAM = ""
      + SAM + LS
      + "4" + TAB + "0" + TAB + "b2" + TAB + "1" + TAB + "255" + TAB + "4=" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ACGT" + TAB + "*" + TAB + "NH:i:2" + TAB + "IH:i:2" + LS
      + "4" + TAB + "0" + TAB + "b2" + TAB + "1" + TAB + "255" + TAB + "4=" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ACGT" + TAB + "*" + TAB + "NH:i:2" + TAB + "IH:i:2" + LS
      + "5" + TAB + "4" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "*" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ACGT" + TAB + "*";

  public void testEndToEnd() throws IOException {
    final File dir = FileUtils.createTempDir("test", "speciesCliDir");
    try {
      final File t = new File(dir, "t");
      ReaderTestUtils.getReaderDNA(">a1\nacgt\n>a2\nacgt\n>b1\nacgt\n>b2\nacgt\n", t, new SdfId(1L)).close();
      final File o = new File(dir, "o");
      final File r = new File(dir, "r");
      FileUtils.stringToFile("a1\ta n a\na2\ta n a\nb1\tb\nb2\tb", r);
      final File s = new File(dir, "s.sam");
      FileUtils.stringToFile(EXTRA_SAM, s);
      final String stats = checkMainInitOk("-o", o.getPath(), "-t", t.getPath(), "-r", r.getPath(), s.getPath(), "--min-confidence", "0.0");
      final String expStats = "Diversity metrics" + LS
          + "       Richness:   2.0000" + LS
          + "        Shannon:  0.68291" + LS
          + "         Pielou:  0.98523" + LS
          + "Inverse Simpson:   1.9600" + LS;
      //System.err.println(stats);
      assertEquals(expStats, stats);
      //System.err.println(stats);
      final String output = FileUtils.fileToString(new File(o, "species.tsv"));
      assertTrue(output.contains("species v2.1"));
      assertTrue(output.contains(SpeciesTask.SPECIES_HEADER));
      //I am unsure of the numbers below
      //System.err.println(output);
      assertTrue(output, output.contains("0.5714 0.1275 1.000 0.4762 0.1063 1.000 2.3 1.500 1.000 8 3.00 Y 1 3 1 species b".replaceAll(" ", "\t")));
      assertTrue(output.contains("0.4286 0.07582 1.000 0.3571 0.06319 1.000 1.6 1.000 1.000 8 2.00 Y 1 2 1 species".replaceAll(" ", "\t") + "\t" + "a n a"));

      final String krona = StringUtils.grepMinusV(FileUtils.fileToString(new File(o, "index.html")), "speciesCliDir");
      mNano.check("speciescli-krona.out", krona);
      final String usage = mCli.usageLog();
      //System.err.println(usage);
      TestUtils.containsAll(usage, "Usage beginning module=species runId=", ", Usage end module=species runId=", " metric=7 success=true");
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  public void testEndToEndFullNames() throws IOException {
    final File dir = FileUtils.createTempDir("test", "speciesCli");
    try {
      final File t = new File(dir, "t");
      ReaderTestUtils.getReaderDNA(">a1 extra1\nacgt\n>a2 stuff2\nacgt\n>b1 before3\nacgt\n>b2 after4\nacgt\n", t, new SdfId(1L)).close();
      final File o = new File(dir, "o");
      final File r = new File(dir, "r");
      FileUtils.stringToFile("a1\ta n a\nb2\tb", r);
      final File s = new File(dir, "s.sam");
      FileUtils.stringToFile(SAM, s);
      final String stats = checkMainInitOk("-o", o.getPath(), "-t", t.getPath(), "-r", r.getPath(), s.getPath(), "--min-confidence", "0.0");
      final String expStats = "Diversity metrics" + LS
          + "       Richness:  4.0000" + LS
          + "        Shannon:  1.3863" + LS
          + "         Pielou:  1.0000" + LS
          + "Inverse Simpson:  4.0000" + LS;
      //System.err.println(stats);
      assertEquals(expStats, stats);
      final String output = FileUtils.fileToString(new File(o, "species.tsv"));
      //System.err.println(output);
      TestUtils.containsAll(output,
          "#TEMPLATE-SDF-ID\t0",
          SpeciesTask.SPECIES_HEADER,
          "0.2500 0.02997 1.000 0.2500 0.02997 1.000 0.88 1.000 1.000 4 1.00 Y 1 2 1 species".replaceAll(" ", "\t") + "\t" + "a n a",
          "0.2500 0.02997 1.000 0.2500 0.02997 1.000 0.88 1.000 1.000 4 1.00 Y 1 4 1 species".replaceAll(" ", "\t") + "\t" + "a2 stuff2",
          "0.2500 0.02997 1.000 0.2500 0.02997 1.000 0.88 1.000 1.000 4 1.00 Y 1 5 1 species".replaceAll(" ", "\t") + "\t" + "b1 before3",
          "0.2500 0.02997 1.000 0.2500 0.02997 1.000 0.88 1.000 1.000 4 1.00 Y 1 3 1 species b".replaceAll(" ", "\t"));

      final String usage =  mCli.usageLog();
      //System.err.println(usage);
      TestUtils.containsAll(usage, "Usage beginning module=species runId=", ", Usage end module=species runId=", " metric=4 success=true");

    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  public void testMyInteger() {
    final SpeciesTask.MyInteger mi = new SpeciesTask.MyInteger();
    assertEquals(1, mi.mValue);
  }

  public void testgetRenameMapFail() throws IOException {
    final File dir = FileUtils.createTempDir("speceiscli", "renamemap");
    final String seq = ">test name" + LS + "acgtacgtgtg" + LS + ">test2 fullname" + LS + "catactgctatgac" + LS;
    final File sdf = ReaderTestUtils.getDNADir(seq, dir);
    try (SequencesReader createDefaultSequencesReader = SequencesReaderFactory.createDefaultSequencesReader(sdf)) {
      String invalidLine = "invalid line with spaces but no tabs";
      final String reference = "#comment line" + LS + "abcd\tlong name with spaces" + LS + invalidLine + LS;
      File refMap = new File(dir, "refmap.txt");
      FileUtils.stringToFile(reference, refMap);
      try {
        SpeciesTask.getRenameMap(refMap, createDefaultSequencesReader);
        fail();
      } catch (final NoTalkbackSlimException ex) {
        assertEquals("The input file " + refMap.getPath() + " for relabel-species-file flag contains invalid entry. line: " + invalidLine, ex.getMessage());
      }
      invalidLine = "\tinvalid";
      refMap = new File(dir, "refmap.txt");
      FileUtils.stringToFile(invalidLine, refMap);
      try {
        SpeciesTask.getRenameMap(refMap, createDefaultSequencesReader);
        fail();
      } catch (final NoTalkbackSlimException ex) {
        assertEquals("The input file " + refMap.getPath() + " for relabel-species-file flag contains invalid entry. line: " + invalidLine, ex.getMessage());
      }
      invalidLine = "invalid\t";
      refMap = new File(dir, "refmap.txt");
      FileUtils.stringToFile(invalidLine, refMap);
      try {
        SpeciesTask.getRenameMap(refMap, createDefaultSequencesReader);
        fail();
      } catch (final NoTalkbackSlimException ex) {
        assertEquals("The input file " + refMap.getPath() + " for relabel-species-file flag contains invalid entry. line: " + invalidLine, ex.getMessage());
      }
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  public void testgetRenameMapWarnings() throws IOException {
    final MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());
    final File dir = FileUtils.createTempDir("speceiscli", "renamemap");
    try {
      final String reference = "#comment line" + LS
          + "abcd\tlong name with spaces" + LS
          + "name2\tname with more         spaces" + LS
          + "name3\tname 3" + LS
          + "name4\tname 4" + LS
          + "name5\tname 5" + LS
          + "name6\tname 6" + LS;
      final String seq = ">test name" + LS + "acgtacgtgtg" + LS + ">test2 fullname" + LS + "catactgctatgac" + LS;
      final File sdf = ReaderTestUtils.getDNADir(seq, dir);
      final File refMap = new File(dir, "refmap.txt");
      FileUtils.stringToFile(reference, refMap);
      try (SequencesReader createDefaultSequencesReader = SequencesReaderFactory.createDefaultSequencesReader(sdf)) {
        final HashMap<String, String> renameMap = SpeciesTask.getRenameMap(refMap, createDefaultSequencesReader);
        assertEquals("long name with spaces", renameMap.get("abcd"));
        assertEquals("name with more         spaces", renameMap.get("name2"));
        assertEquals("test name", renameMap.get("test"));
        assertEquals("test2 fullname", renameMap.get("test2"));
        //System.err.println(mps.toString());
        TestUtils.containsAll(mps.toString(), "Template SDF does not contain any sequence named \"abcd\" specified in the " + refMap.getPath() + " file",
          "Template SDF does not contain any sequence named \"name2\" specified in the " + refMap.getPath() + " file",
          "Template SDF does not contain any sequence named \"name4\" specified in the " + refMap.getPath() + " file",
          "Template SDF does not contain any sequence named \"name5\" specified in the " + refMap.getPath() + " file",
          "Template SDF does not contain any sequence named \"name5\" specified in the " + refMap.getPath() + " file",
          "Subsequent warnings of this type will not be shown.",
          "There were 6 names not present in the template SDF.");
      }

    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  public void testSortWarning() throws IOException {
    final File dir = FileUtils.createTempDir("test", "speciesCli");
    try {
      final File t = new File(dir, "t");
      ReaderTestUtils.getReaderDNA(">a1\nacgt\n>a2\nacgt\n>b1\nacgt\n>b2\nacgt\n", t, new SdfId(1L)).close();
      final File o = new File(dir, "o");
      final File r = new File(dir, "r");
      FileUtils.stringToFile("a1\ta n a\na2\ta n a\nb1\tb\nb2\tb", r);
      final File s = new File(dir, "s.sam");
      FileUtils.stringToFile(SharedSamConstants.SAM_UNSORTED, s);
      final String[] args = {
          "species", "-o", o.getPath(), "-t", t.getPath(), "-r", r.getPath(), s.getPath()
      };
      final MemoryPrintStream mps = new MemoryPrintStream();
      try {
        new Slim().intMain(args, mps.outputStream(), mps.printStream());
      } finally {
        mps.close();
      }
      //System.err.println(mps.toString());
      assertTrue(mps.toString().contains("is not sorted in coordinate order."));
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }
}

