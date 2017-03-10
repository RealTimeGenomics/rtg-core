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
package com.rtg.ngs;

import java.io.IOException;

import com.rtg.ngs.blocking.MapQScoringReadBlocker;
import com.rtg.ngs.blocking.MapQScoringReadBlockerSynch;
import com.rtg.reader.NamesInterface;
import com.rtg.reader.PrereadType;
import com.rtg.reader.SequencesReader;
import com.rtg.reference.Sex;
import com.rtg.sam.SamUtils;
import com.rtg.util.Constants;
import com.rtg.util.Environment;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.NoTalkbackSlimException;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMSequenceRecord;

/**
 * Resources that can be shared among multiple instances of
 * <code>SamPairedAlignmentWriter</code>
 *
 */
public class SharedResources {

  /**
   * Create shared resources.
   *
   * @param params parameters shared by multiple instances of
   *          <code>SamPairedAlignmentWriter</code>
   * @return the resources object
   * @throws IOException if an exception occurs during IO
   */
  public static SharedResources generateSharedResources(final NgsParams params) throws IOException {
    final SequencesReader srFirst = params.buildFirstParams().reader().copy();
    final SequencesReader srSecond = params.buildSecondParams() != null ? params.buildSecondParams().reader().copy() : null;
    final SequencesReader srTemplate = params.searchParams().reader().copy();

    assert params.outputParams() != null;
    final MapQScoringReadBlocker blocker = new MapQScoringReadBlockerSynch((int) srFirst.numberSequences(), params.outputParams().maxTopResults());
    final NamesInterface templateNames = params.searchParams().reader().names();
    final SAMFileHeader header = createHeader(srTemplate, templateNames, srFirst, params.outputParams().readGroup(), params.sex(), true);
    final SAMFileHeader header2 = createHeader(srTemplate, templateNames, srFirst, params.outputParams().readGroup(), params.sex(), false);
    final SharedResources sr = new SharedResources(srFirst, srSecond, srTemplate, blocker, templateNames, header, header2);
    if (params.useTopRandom()) {
      PairedTopRandomImplementation pairedEndTopRandom = null;
      final SingleEndTopRandomImplementation singleEndTopRandom;
      if (params.numberThreads() > 1) {
        if (params.paired()) {
          pairedEndTopRandom = new PairedTopRandomImplementationSync((int) srFirst.numberSequences());
          singleEndTopRandom = new SingleEndTopRandomImplementationSync((int) srFirst.numberSequences() * 2);
        } else {
          singleEndTopRandom = new SingleEndTopRandomImplementationSync((int) srFirst.numberSequences());
        }
      } else {
        if (params.paired()) {
          pairedEndTopRandom = new PairedTopRandomImplementation((int) srFirst.numberSequences());
          singleEndTopRandom = new SingleEndTopRandomImplementation((int) srFirst.numberSequences() * 2);
        } else {
          singleEndTopRandom = new SingleEndTopRandomImplementation((int) srFirst.numberSequences());
        }
      }
      if (params.paired()) {
        pairedEndTopRandom.initialize();
      }
      sr.setTopRandom(pairedEndTopRandom, singleEndTopRandom);
    }
    return sr;
  }


  private final SequencesReader mFirstReader;
  private final SequencesReader mSecondReader;
  private final SequencesReader mTemplateReader;
  private final MapQScoringReadBlocker mBlocker;
  private final NamesInterface mTemplateNames;
  private final SAMFileHeader mFileHeader;
  private PairedTopRandomImplementation mPairedEndTopRandom = null;
  private SingleEndTopRandomImplementation mSingleEndTopRandom = null;

  SharedResources(SequencesReader first, SequencesReader second, SequencesReader template,
                  MapQScoringReadBlocker blocker, NamesInterface templateNames, SAMFileHeader header, SAMFileHeader headerNoDict) {
    mFirstReader = first;
    mSecondReader = second;
    mTemplateReader = template;
    validateReader(mFirstReader);
    validateReader(mSecondReader);
    mBlocker = blocker;
    mTemplateNames = templateNames;
    mFileHeader = header;
  }

  protected static SAMFileHeader createHeader(SequencesReader templateReader, NamesInterface templateNames, SequencesReader leftReader, SAMReadGroupRecord readGroupRecord, Sex sex, boolean includeDictionary) throws IOException {
    final SAMFileHeader header = new SAMFileHeader();
    header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
    final SAMProgramRecord pg = new SAMProgramRecord(Constants.APPLICATION_NAME);
    pg.setProgramVersion(Environment.getVersion());
    header.addComment(SamUtils.RUN_ID_ATTRIBUTE + CommandLine.getRunId());
    if (templateReader.getSdfId().available()) {
      header.addComment(SamUtils.TEMPLATE_SDF_ATTRIBUTE + templateReader.getSdfId().toString());
    }
    if (leftReader.getSdfId().available()) {
      header.addComment(SamUtils.READ_SDF_ATTRIBUTE + leftReader.getSdfId().toString());
    }
    if (sex != null) {
      header.addComment(SamUtils.GENDER_ATTRIBUTE + sex.toString());
    }
    if (CommandLine.getCommandLine() != null) {
      pg.setCommandLine(CommandLine.getCommandLine());
    }
    SamUtils.addProgramRecord(header, pg);
    if (includeDictionary) {
      final int[] lengths = templateReader.sequenceLengths(0, templateReader.numberSequences());
      for (int i = 0; i < lengths.length; ++i) {
        final SAMSequenceRecord record = new SAMSequenceRecord(templateNames.name(i), lengths[i]);
        header.addSequence(record);
      }
    }
    if (readGroupRecord != null) {
      header.addReadGroup(readGroupRecord);
    }
    return header;
  }

  void setTopRandom(PairedTopRandomImplementation pairedEndTopRandom, SingleEndTopRandomImplementation singleEndTopRandom) {
    mPairedEndTopRandom = pairedEndTopRandom;
    mSingleEndTopRandom = singleEndTopRandom;
  }

  private void validateReader(final SequencesReader reader) {
    if (reader != null && reader.getPrereadType() == PrereadType.CG && reader.minLength() != reader.maxLength()) {
      throw new NoTalkbackSlimException(ErrorType.CG_LENGTH_ERROR, reader.path() == null ? "input tsv file" : reader.path().getAbsolutePath());
    }
  }

  /**
   * Test if this run is paired-end.
   *
   * @return true for paired-end
   */
  public boolean isPairedEnd() {
    return mSecondReader != null;
  }

  public SAMFileHeader getHeader() {
    return mFileHeader;
  }

  public MapQScoringReadBlocker getBlocker() {
    return mBlocker;
  }

  public PairedTopRandomImplementation getPairedEndTopRandom() {
    return mPairedEndTopRandom;
  }
  public SingleEndTopRandomImplementation getSingleEndTopRandom() {
    return mSingleEndTopRandom;
  }

  /**
   * Return the names of all template sequences.
   * @return template names
   */
  public NamesInterface names() {
    return mTemplateNames;
  }

  /**
   * @return a copy of the first reader
   */
  public SequencesReader firstReaderCopy() {
    return mFirstReader == null ? null : mFirstReader.copy();
  }

  /**
   * @return a copy of the second reader
   */
  public SequencesReader secondReaderCopy() {
    return mSecondReader == null ? null : mSecondReader.copy();
  }

  /**
   * @return a copy of the template reader
   */
  public SequencesReader templateReaderCopy() {
    return mTemplateReader == null ? null : mTemplateReader.copy();
  }

  /**
   * Closes all the resources
   * @throws IOException if error
   */
  public void close() throws IOException {
    if (mFirstReader != null) {
      mFirstReader.close();
    }
    if (mSecondReader != null) {
      mSecondReader.close();
    }
    if (mTemplateReader != null) {
      mTemplateReader.close();
    }
    if (mBlocker != null) {
      mBlocker.close();
    }
  }
}
