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

package com.rtg.assembler;

import java.io.File;
import java.io.IOException;
import java.util.Collections;

import com.rtg.assembler.graph.MutableGraph;
import com.rtg.assembler.graph.implementation.GraphKmerAttribute;
import com.rtg.assembler.graph.io.GraphWriter;
import com.rtg.launcher.AbstractParamsCliTest;
import com.rtg.launcher.ParamsCli;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SdfId;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.cli.CFlags;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.store.StoreDirProxy;
import com.rtg.util.test.FileHelper;

/**
 */
public class AssembleCliTest extends AbstractParamsCliTest<AssembleParams> {
  @Override
  protected ParamsCli<AssembleParams> getParamsCli() {
    return new AssembleCli();
  }

  public void testFlags() throws IOException {
    final MutableGraph builtGraph = GraphMapCliTest.makeGraph(29
        , new String[]{"ACGT", "GGGG", "TTAA"}
        , new long[][]{{1, 2}, {3, 2}});
    final File tmpDir = FileHelper.createTempDirectory();
    try {
      final File graph = new File(tmpDir, "input");
      final File reads = ReaderTestUtils.getDNADir(">a" + StringUtils.LS + "ACGTACGTACGTACGTACGTACGTACGTAC" + StringUtils.LS);
      assertTrue(graph.mkdir());
      try {
        GraphWriter.write(builtGraph, new StoreDirProxy(graph), "monkey", Collections.emptySet());

        try {
          final File fileList = new File(tmpDir, "fileList");
          FileUtils.stringToFile(reads.toString() + StringUtils.LS, fileList);

          checkHandleFlagsErr();
          final File output = new File(tmpDir, "bar");
          String err = checkHandleFlagsErr("-o", output.toString(), reads.toString());
          TestUtils.containsAllUnwrapped(err, "You must provide a value for -k INT");
          err = checkHandleFlagsErr(reads.toString(), "-g", graph.toString(), "-k", "30");
          TestUtils.containsAllUnwrapped(err, "You must provide a value for -o DIR");

          err = checkHandleFlagsErr("-o", output.toString(), reads.toString(), "-g", graph.toString(), "-w", "-1", "-k", "30", "--consensus-reads", "10");
          TestUtils.containsAllUnwrapped(err, "--word", "must be at least 1");
          err = checkHandleFlagsErr("-o", output.toString(), reads.toString(), "-g", graph.toString(), "-s", "-1", "-k", "30", "--consensus-reads", "10");
          TestUtils.containsAllUnwrapped(err, "--step", "must be at least 1");

          err = checkHandleFlagsErr("-o", output.toString(), reads.toString(), "-g", graph.toString(), "-k", "-30", "--consensus-reads", "10");
          TestUtils.containsAllUnwrapped(err, "--kmer-size", "must be at least 1");

          err = checkHandleFlagsErr("-o", output.toString(), reads.toString(), "-g", graph.toString(), "-k", "30", "-a", "-1", "--consensus-reads", "10");
          TestUtils.containsAllUnwrapped(err, "The specified flag \"--mismatches\" has invalid value \"-1\"");

          err = checkHandleFlagsErr("-o", output.toString(), reads.toString(), "-g", graph.toString(), "-k", "30", "-m", "15", "-M", "10", "--consensus-reads", "10");
          TestUtils.containsAllUnwrapped(err, "--max-insert should be larger than --min-insert");

          checkHandleFlags("-o", output.toString(), reads.toString(), "-g", graph.toString(), "-k", "30", "-a", "10%");
          checkHandleFlags("-o", output.toString(), reads.toString(), "-g", graph.toString(), "-k", "30");
          FileHelper.deleteAll(output);
          checkMainInitOk("-o", output.toString(), "-f", reads.toString(), "-k", "3", "--consensus-reads", "10");
          FileHelper.deleteAll(output);
          checkMainInitOk("-o", output.toString(), "-F", fileList.toString(), "-k", "3", "--consensus-reads", "10");
          FileHelper.deleteAll(output);
          checkMainInitOk("-o", output.toString(), "-j", reads.toString(), "-k", "3", "--consensus-reads", "10");
          FileHelper.deleteAll(output);
          checkMainInitOk("-o", output.toString(), "-J", fileList.toString(), "-k", "3", "--consensus-reads", "10");
          FileHelper.deleteAll(output);
          checkMainInitOk("-o", output.toString(), "-J", fileList.toString(), "-g", graph.toString(), "-k", "30", "--consensus-reads", "10");
          FileHelper.deleteAll(output);
          checkMainInitOk("-o", output.toString(), "-I", fileList.toString(), "-g", graph.toString(), "-k", "30", "--consensus-reads", "10");

          final File inputGraph =  new File(tmpDir, "graph");
          assertTrue(inputGraph.mkdir());
          GraphWriter.write(new GraphKmerAttribute(30), new StoreDirProxy(inputGraph), "foo" , Collections.emptySet());
          final MemoryPrintStream mps = new MemoryPrintStream();
          CFlags flags =  new CFlags("foo", mps.printStream(), mps.printStream());
          AssembleCli.initLocalFlags(flags);
          assertFalse(mps.toString(), flags.setFlags("-o", output.toString(), "-g", inputGraph.toString(), "-I", fileList.toString(), "-w", "29", "-s", "12", "-k", "30", "--consensus-reads", "10"));
          FileHelper.deleteAll(output);
          assertTrue(mps.toString(), flags.setFlags("-o", output.toString(), "-g", inputGraph.toString(), "-I", fileList.toString(), "-w", "29", "-s", "12", "-k", "30", "--consensus-reads", "10"));
          AssembleParams params = AssembleCli.makeParamsLocal(flags);
          assertEquals(29, params.wordSize());
          assertEquals(reads, params.reads().get(0));
          assertEquals(12, params.stepSize());
          mps.reset();
          assertTrue(mps.toString(), flags.setFlags("-o", output.toString(), "-I", fileList.toString(), "-g", inputGraph.toString(), "-k", "30", "--consensus-reads", "10", "-p", "13", "-r", "14"));
          params = AssembleCli.makeParamsLocal(flags);
          assertEquals(18, params.wordSize());
          assertEquals(18, params.stepSize());
          assertEquals(10, params.consensusThreshold());
          assertEquals(13, params.minPathReads());
          assertEquals(14, params.minReadCount());
          mps.reset();
          flags =  new CFlags("foo", mps.printStream(), mps.printStream());
          AssembleCli.initLocalFlags(flags);
          assertTrue(mps.toString(), flags.setFlags("-o", output.toString(), "-I", fileList.toString(), "-g", inputGraph.toString(), "-k", "30", "--consensus-reads", "10"));
          params = AssembleCli.makeParamsLocal(flags);
          assertEquals(-1, params.minPathReads());
          assertEquals(-1, params.minReadCount());
        } finally {
          FileHelper.deleteAll(graph);
        }
      } finally {
        FileHelper.deleteAll(reads);
      }
    } finally {
      FileHelper.deleteAll(tmpDir);
    }
  }

  public void testInitParams() {
    checkHelp("assemble [OPTION]... -k INT -o DIR SDF+",
        "-k INT -o DIR -I FILE",
        "-k INT -o DIR -f SDF",
        "-k INT -o DIR -F FILE",
        "--consensus-reads",
        "-g, --graph=DIR ", "graph of the assembly to map against",
        "-w, --word=INT", "word size",
        "-s, --step=INT", "step size",
        "-k, --kmer-size=INT", "kmer length to build graph nodes from",
        "-a, --mismatches=INT", "number of bases that may mismatch in an alignment",
        "-m, --min-insert=INT", "minimum insert size between fragments",
        "-M, --max-insert=INT", "maximum insert size between fragments",
        "-f, --454", "SDF containing 454 reads",
        "-F, --input-list-454", "file containing",
        "-j, --mate-pair", "SDF containing mate pair reads",
        "-J, --input-list-mate-pair",
        "-c, --minimum-kmer-frequency", "set minimum kmer frequency to retain, or -1 for automatic threshold",
        "-r, --read-count", "minimum reads required to overlap a link",
        "--preserve-bubbles", "avoid merging bubbles where the ratio of kmers on the branches is below this",
        "-p, --min-path", "minimum reads required in long paths"
    );
  }

  // A bug arising because the --end-read is much larger than actual number of reads (1)
  // Manifests at a deep exception in AsyncReadSource
  public void testBug1527() throws IOException {
    final File tmpDir = FileHelper.createTempDirectory();
    try {
      final StringBuilder sb = new StringBuilder();
      for (int k = 0; k < 750; ++k) {
        sb.append('A');
      }
      final File reads = ReaderTestUtils.getDNADir(">a" + StringUtils.LS + sb.toString() + StringUtils.LS);
      try {
        final File fileList = new File(tmpDir, "fileList");
        FileUtils.stringToFile(reads.toString() + StringUtils.LS, fileList);
        final File output = new File(tmpDir, "bar");
        final String err = checkMainInitWarn("-o", output.toString(), "--kmer-size", "750", "--end-read", "2147483647", reads.toString());
        assertEquals(err, "The end sequence id \"2147483647\" is out of range, it must be from \"1\" to \"1\". Defaulting end to \"1\"" + StringUtils.LS);
      } finally {
        FileHelper.deleteAll(reads);
      }
    } finally {
      FileHelper.deleteAll(tmpDir);
    }
  }

  public void testBug1527PairedEnd() throws IOException {
    final File tmpDir = FileHelper.createTempDirectory();
    try {
      final StringBuilder sb = new StringBuilder();
      for (int k = 0; k < 750; ++k) {
        sb.append('A');
      }
      final String fasta = ">a" + StringUtils.LS + sb.toString() + StringUtils.LS;
      final File reads = new File(tmpDir, "sdf");
      ReaderTestUtils.createPairedReaderDNA(fasta, fasta, reads, new SdfId());
      try {
        final File fileList = new File(tmpDir, "fileList");
        FileUtils.stringToFile(reads.toString() + StringUtils.LS, fileList);
        final File output = new File(tmpDir, "bar");
        final String err = checkMainInitWarn("-o", output.toString(), "--kmer-size", "750", "--end-read", "2147483647", reads.toString());
        assertEquals(err, "The end sequence id \"2147483647\" is out of range, it must be from \"1\" to \"1\". Defaulting end to \"1\"" + StringUtils.LS);
      } finally {
        FileHelper.deleteAll(reads);
      }
    } finally {
      FileHelper.deleteAll(tmpDir);
    }
  }
}
