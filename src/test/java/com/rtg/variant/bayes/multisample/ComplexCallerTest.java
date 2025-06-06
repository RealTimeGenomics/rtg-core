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

package com.rtg.variant.bayes.multisample;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import com.rtg.launcher.MockReaderParams;
import com.rtg.mode.DnaUtils;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.relation.RelationshipsFileParser;
import com.rtg.sam.CircularBufferMultifileSinglePassReaderWindow;
import com.rtg.sam.CircularBufferMultifileSinglePassReaderWindowTest;
import com.rtg.sam.ReaderWindow;
import com.rtg.sam.RecordIterator;
import com.rtg.sam.SamFilterParams;
import com.rtg.sam.SamRegionRestriction;
import com.rtg.sam.SamUtils;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.Populator;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.intervals.RegionRestriction;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;
import com.rtg.variant.DefaultMachineErrorChooser;
import com.rtg.variant.GenomePriorParamsBuilder;
import com.rtg.variant.StaticThreshold;
import com.rtg.variant.Variant;
import com.rtg.variant.VariantAlignmentRecord;
import com.rtg.variant.VariantAlignmentRecordPopulator;
import com.rtg.variant.VariantOutputLevel;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantParamsBuilder;
import com.rtg.variant.bayes.multisample.ComplexRegion.RegionType;
import com.rtg.variant.bayes.multisample.cancer.SomaticCallerConfiguration;
import com.rtg.variant.bayes.multisample.singleton.SingletonCallerConfiguration;
import com.rtg.variant.format.VariantOutputVcfFormatter;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import junit.framework.TestCase;

/**
 */
public class ComplexCallerTest extends TestCase {

  private SAMFileHeader mHeader = null;
  private static final String TEMPLATE_NAME = "TEST";
  private static final String TEMPLATE = "ACGTACGTACGTACGTACGTACGTACGT";

  @Override
  protected void setUp() {
    Diagnostic.setLogStream();
    mHeader = new SAMFileHeader();
    mHeader.addSequence(new SAMSequenceRecord("foo", 0));
  }
  @Override
  protected void tearDown() {
    mHeader = null;
  }

  private static void addSample(final SAMFileHeader header, final String sampleName) {
    final SAMReadGroupRecord rg = new SAMReadGroupRecord(sampleName);
    rg.setSample(sampleName);
    header.addReadGroup(rg);
  }

  public static SAMFileHeader makeHeaderWithSamples(final String... sampleNames) {
    final SAMFileHeader uber = new SAMFileHeader();
    for (final String sampleName : sampleNames) {
      addSample(uber, sampleName);
    }
    return uber;
  }

  public void testComplexCaller() throws Exception {
    final VariantParamsBuilder builder = new VariantParamsBuilder();
    builder.genomePriors(new GenomePriorParamsBuilder().create());
    builder.machineErrorName("default");
    builder.maxCoverageFilter(new StaticThreshold(100));
    builder.callLevel(VariantOutputLevel.ALL);
    builder.genomeRelationships(RelationshipsFileParser.load(new BufferedReader(new StringReader("original-derived normal cancer contamination=0.5"))));
    final SAMFileHeader header = makeHeaderWithSamples("normal", "cancer");
    final VariantParams params = builder
      .uberHeader(header)
      .genome(new MockReaderParams(ReaderTestUtils.getReaderDnaMemory(ReaderTestUtils.SEQ_DNA_SIMPLE)))
      .create();
    final ArrayList<Variant> chunk = new ArrayList<>();
    chunk.add(TestUtils.createVariant(4));
    chunk.add(TestUtils.createVariant(5));
    chunk.add(TestUtils.createVariant(6));
    chunk.add(TestUtils.createVariant(10, true));
    chunk.add(TestUtils.createVariant(15));
    final MockReaderWindow trib = new MockReaderWindow();
    final Complexities regions = new Complexities(chunk, "foo", 0, 50, 3, 15, ComplexitiesTest.template(30), true, null);
    Complexities.fixDangling(regions, null);
    Complexities.fixDangling(null, regions);
    final AbstractJointCallerConfiguration config = new SomaticCallerConfiguration.Configurator().getConfig(params, null);
    final ComplexCaller caller = new ComplexCaller(params, config);
    final List<Variant> list = caller.makeComplexCalls(regions, trib, DnaUtils.encodeString(TEMPLATE), TEMPLATE_NAME);
    assertEquals(list.toString(), 2, list.size());

    assertEquals(4, list.get(0).getLocus().getStart());
    assertEquals(10, list.get(1).getLocus().getStart());
  }

  private class MockReaderWindow implements ReaderWindow<VariantAlignmentRecord> {

    @Override
    public Iterator<VariantAlignmentRecord> recordsOverlap(int start, int end) {
      return new Iterator<VariantAlignmentRecord>() {

        private int mCount = 0;

        @Override
        public boolean hasNext() {
          return mCount++ <= 20;
        }

        @Override
        public VariantAlignmentRecord next() {
          final SAMRecord rec = new SAMRecord(mHeader);
          rec.setAlignmentStart(mCount);
          rec.setReadName("" + mCount);
          rec.setCigarString("3M1I4M");
          rec.setReadString("CCGTCCGT");
          rec.setReferenceIndex(0);
          rec.setMappingQuality(10);
          return new VariantAlignmentRecord(rec);
        }

        @Override
        public void remove() {
        }
      };
    }

    @Override
    public void advanceBuffer(int end) {
      // do nothing
    }

    @Override
    public void flush(int start, int end) {
      // do nothing
    }
    @Override
    public int flushedTo() {
      return 0;
    }
  }

   private class MockComplexities extends Complexities {

    MockComplexities() {
      super(new ArrayList<>(), "foo", 0, 2000, 5, 5, ComplexitiesTest.template(30), true, null);
    }

    @Override
    public boolean isFixed() {
      return true;
    }

    @Override
    public Iterator<ComplexRegion> iterator() {
      return new Iterator<ComplexRegion>() {

        int mCount = 0;

        @Override
        public boolean hasNext() {
          return mCount == 0;
        }

        @Override
        public ComplexRegion next() {
          if (mCount == 0) {
            ++mCount;
            return new ComplexRegion("foo", 47, 53, RegionType.INTERESTING);
          }
          return null;
        }

        @Override
        public void remove() {
        }
      };
    }
  }

  public void testComplexCallReadOrders() throws IOException, InvalidParamsException {
    try (final TestDirectory dir = new TestDirectory("complexcall-readorder")) {
      final File in = new File(dir, "complexreadorder.sam.gz");
      FileHelper.resourceToFile("com/rtg/variant/resources/complexreadorder.sam.gz", in);
      final File index = new File(dir, "complexreadorder.sam.gz.tbi");
      FileHelper.resourceToFile("com/rtg/variant/resources/complexreadorder.sam.gz.tbi", index);
      check(in);
    }
  }

  private class MyCB extends CircularBufferMultifileSinglePassReaderWindow<VariantAlignmentRecord> {

    private List<SAMRecord> readFile(File in) throws IOException {
      final List<SAMRecord> list = new ArrayList<>();
      try (SamReader r = SamUtils.makeSamReader(in)) {
        for (Object aR : r) {
          list.add((SAMRecord) aR);
        }
      }
      return list;
    }

    private final List<SAMRecord> mList;

    MyCB(RecordIterator<VariantAlignmentRecord> recordIt, List<File> samFiles, RegionRestriction restriction, Populator<VariantAlignmentRecord> pop, SAMFileHeader header) throws IOException {
      super(recordIt, pop, header.getSequenceIndex("chr21"), restriction.getStart(), Integer.MAX_VALUE);
      assert samFiles.size() == 1;
      mList = readFile(samFiles.get(0));
      //System.err.println("size= " + mList.size());
    }

    @Override
    public Iterator<VariantAlignmentRecord> recordsOverlap(int start, int end) {
      return new Iterator<VariantAlignmentRecord>() {

        private int mStart = mList.size() - 1;

        @Override
        public boolean hasNext() {
          return mStart >= 0;
        }

        @Override
        public VariantAlignmentRecord next() {
          return new VariantAlignmentRecord(mList.get(mStart--));
        }

        @Override
        public void remove() {
          throw new UnsupportedOperationException();
        }
      };
    }
  }

  /**
   * @throws InvalidParamsException
   * @throws IOException
   */
  private void check(File in) throws InvalidParamsException, IOException {
    final VariantParamsBuilder builder = new VariantParamsBuilder();
    builder.genomePriors(new GenomePriorParamsBuilder().create());
    builder.machineErrorName("illumina");
    builder.maxCoverageFilter(new StaticThreshold(300));
    final String template = "AGCATTTTTGAAATTCTCTTTTTGTAATATCTGCAAGTAGACATTTGGAGTACTTTGAGGCCTATTGTGGAAAAGGAAATATCTTCACAGAAAAACTAGATA";
    final SequencesReader reader = ReaderTestUtils.getReaderDnaMemory(">chr21\n" + template);
    final SamRegionRestriction restriction = new SamRegionRestriction("chr21", 0, template.length());
    builder.filterParams(SamFilterParams.builder().findAndRemoveDuplicates(false).restriction(restriction).create());

    final List<File> list = new ArrayList<>();
    list.add(in);
    //System.err.println(list.size() + " " + in.getPath());
    final SAMFileHeader header = SamUtils.getUberHeader(list);
    final VariantParams params = builder
      .uberHeader(header)
      .genome(new MockReaderParams(reader))
      .create();
    final AbstractJointCallerConfiguration config = new SingletonCallerConfiguration.Configurator().getConfig(params, null);

    final VariantAlignmentRecordPopulator pop = new VariantAlignmentRecordPopulator(new DefaultMachineErrorChooser(), 0);
    final RecordIterator<VariantAlignmentRecord> it = CircularBufferMultifileSinglePassReaderWindowTest.defaultIterator(list, params.filterParams(), pop);
    final CircularBufferMultifileSinglePassReaderWindow<VariantAlignmentRecord> trib =
       new CircularBufferMultifileSinglePassReaderWindow<>(
           it,
            pop,
            header.getSequenceIndex("chr21"), restriction.getStart(), Integer.MAX_VALUE
            );
    final ComplexCaller caller = new ComplexCaller(params, config);
    final List<Variant> res = caller.makeComplexCalls(new MockComplexities(), trib, DnaUtils.encodeString(template), "chr21");
    final RecordIterator<VariantAlignmentRecord> it2 = CircularBufferMultifileSinglePassReaderWindowTest.defaultIterator(list, params.filterParams(), pop);
    final CircularBufferMultifileSinglePassReaderWindow<VariantAlignmentRecord> trib2 =
        new MyCB(it2, list,
             restriction,
             pop,
            header
        );
     final List<Variant> res2 = caller.makeComplexCalls(new MockComplexities(), trib2, DnaUtils.encodeString(template), "chr21");
     assertEquals(res.size(), res2.size());
     final VariantOutputVcfFormatter format = new VariantOutputVcfFormatter(params, "test");
     assertEquals(format.formatCall(res.get(0)), format.formatCall(res2.get(0)));
    trib.close();
    it.close();
    trib2.close();
    it2.close();
  }
}
