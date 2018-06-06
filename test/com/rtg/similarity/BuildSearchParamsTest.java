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
package com.rtg.similarity;

import static com.rtg.util.StringUtils.LS;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.rtg.index.params.CountParams;
import com.rtg.index.params.CreateParams;
import com.rtg.launcher.BuildParams;
import com.rtg.launcher.HashingRegion;
import com.rtg.launcher.SequenceParams;
import com.rtg.mode.ProgramMode;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.Pair;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 */
public class BuildSearchParamsTest extends TestCase {

  protected File mDir;

  @Override
  public void setUp() throws IOException {
    Diagnostic.setLogStream();
    mDir = FileHelper.createTempDirectory();
  }

  @Override
  public void tearDown() {
    assertTrue(FileHelper.deleteAll(mDir));
    mDir = null;
  }

  BuildSearchParams getParams(final ProgramMode mode, final BuildParams build, final CountParams count) throws IOException {
    return BuildSearchParams.builder().mode(mode).build(build).index(CreateParams.fromBuildParams(build).create())
                                      .count(count).create();
  }

  BuildSearchParams getParams(final ProgramMode mode, final SequenceParams sequence, final BuildParams build, final CountParams count, boolean emptySequenceList) {
    final ArrayList<Pair<String, List<SequenceParams>>> sequences = new ArrayList<>();
    if (!emptySequenceList) {
      final Pair<String, List<SequenceParams>> pair = new Pair<>("test", new ArrayList<>());
      pair.getB().add(sequence);
      sequences.add(pair);
    }
    return BuildSearchParams.builder().mode(mode).sequences(sequences).name("SequencesBuildParams")
                                      .uniqueWords(true).count(count).build(build).create();
  }

  public void testEquals() throws IOException, ClassNotFoundException {
    Diagnostic.setLogStream();
    final ProgramMode pma = ProgramMode.SLIMN;
    final ProgramMode pmb = ProgramMode.TSLIMX;

    final File subjectDir = ReaderTestUtils.getDNADir(mDir);
    final SequenceParams subjectaa = SequenceParams.builder().directory(subjectDir).mode(pma.subjectMode()).create();
    final BuildParams buildaa = BuildParams.builder().windowSize(4).stepSize(1).sequences(subjectaa).create();
    final SequenceParams subjectab = SequenceParams.builder().directory(subjectDir).mode(pmb.subjectMode()).create();
    final BuildParams buildab = BuildParams.builder().windowSize(4).stepSize(2).sequences(subjectab).create();

    final SequenceParams subjectb = SequenceParams.builder().directory(subjectDir).mode(pmb.subjectMode()).create();
    final BuildParams buildb = BuildParams.builder().windowSize(4).stepSize(1).sequences(subjectb).create();

    final File hitDir = FileHelper.createTempDirectory();
    try {
      final CountParams count = new CountParams(hitDir, 5, 10, false);

      final BuildSearchParams a = getParams(pma, buildaa, count);
      final BuildSearchParams a1 = getParams(pma, buildaa, count);
      final BuildSearchParams b = getParams(pmb, buildb, count);
      final BuildSearchParams c = getParams(pmb, buildab, count);
      TestUtils.equalsHashTest(new BuildSearchParams[][] {{a, a1}, {b}, {c}});
      a.close();
      b.close();
      c.close();
    } finally {
      assertTrue(FileHelper.deleteAll(hitDir));
    }
  }

  public void testSequencesList() throws Exception {
    Diagnostic.setLogStream();
    final ProgramMode pm = ProgramMode.SLIMN;
    final SequenceParams dummySubjectParams = SequenceParams.builder().region(new HashingRegion(0, 1)).mode(pm.subjectMode()).create();
    final File hitDir = FileHelper.createTempDirectory();
    try {
      final CountParams count = new CountParams(hitDir, 5, 10, false);
      final File subjectDir = ReaderTestUtils.getDNADir(mDir);
      final SequenceParams subject = SequenceParams.builder().directory(subjectDir).mode(pm.subjectMode()).create();
      final BuildParams build = BuildParams.builder().windowSize(4).stepSize(1).sequences(dummySubjectParams).create();


      try (BuildSearchParams bsp = getParams(pm, subject, build, count, false)) {
        bsp.integrity();
        assertEquals(4, bsp.bufferLength());
        assertTrue(bsp.sequences() instanceof List);
      }

      try (BuildSearchParams bsp = getParams(pm, subject, build, count, true)) {
        bsp.integrity();
        assertEquals(0, bsp.bufferLength());
        assertTrue(bsp.sequences() instanceof List);
      }
    } finally {
      assertTrue(FileHelper.deleteAll(hitDir));
    }
  }

  public void test() throws Exception {
    Diagnostic.setLogStream();
    final ProgramMode pm = ProgramMode.SLIMN;
    final File subjectDir = ReaderTestUtils.getDNADir(mDir);
    final SequenceParams subject = SequenceParams.builder().directory(subjectDir).mode(pm.subjectMode()).create();
    final BuildParams build = BuildParams.builder().windowSize(4).stepSize(1).sequences(subject).create();

    final File hitDir = FileHelper.createTempDirectory();
    try {
      final CountParams count = new CountParams(hitDir, 5, 10, false);

      try (BuildSearchParams bsp = getParams(pm, build, count)) {
        bsp.integrity();

        assertEquals(pm, bsp.mode());
        assertEquals(build.toString(), bsp.build().toString());
        assertEquals(4, bsp.bufferLength());
        assertFalse(bsp.uniqueWords());
        assertEquals(hitDir, bsp.directory());
        assertEquals(new File(hitDir, "bob"), bsp.file("bob"));
        assertEquals(""
            + "BuildSearchParams mode=SLIMN" + LS
            + ".. hits={ CountParams directory=" + bsp.countParams().directory()
            + " topN=5 minHits=10 max. file size=1000000000} " + LS
            + ".. index={ size=1 hash bits=8 initial pointer bits=2 value bits=20}" + LS
            + ".. build={ seq={SequenceParams mode=UNIDIRECTIONAL region=[(0:-1), (1:-1)] directory="
            + build.directory().toString()
            + "}  window=4 step=1}" + LS
          , bsp.toString()
        );
      }
    } finally {
      assertTrue(FileHelper.deleteAll(hitDir));
    }
  }

}

