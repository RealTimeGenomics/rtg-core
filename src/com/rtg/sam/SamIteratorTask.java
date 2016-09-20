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

import java.io.IOException;
import java.io.OutputStream;
import java.util.Map;

import com.rtg.launcher.ParamsTask;
import com.rtg.launcher.Statistics;
import com.rtg.reader.ReaderUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.util.MathUtils;
import com.rtg.util.SingletonPopulatorFactory;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.diagnostic.WarningType;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;

/**
 * Iterator over sorted SAM records. This is mostly deprecated in favour of using the tribble method.
 * @param <P> type of underlying implementation.
 * @param <S> type of statistics object.
 */
public abstract class SamIteratorTask<P extends SingleMappedParams, S extends Statistics> extends ParamsTask<P, S> {

  protected int mPreviousStart = -1;
  protected int mTemplateLength = -2;
  protected String mTemplateName = null;
  private long mNumberTemplateSequences;
  private long mProgressNT;
  private long mTotalTemplateLength;
  static final long PROGRESS_NT = 1000000;
  static final long PROGRESS_NT_DIV = 100;
  protected final SequencesReader mGenomeSequences;
  protected final Map<String, Long> mTemplateNameMap;
  protected final SamFilterParams mFilterParams;

  protected SamIteratorTask(final P params, final OutputStream defaultOutput, final S stats, SamFilterParams filteringParams) {
    super(params, defaultOutput, stats, null);
    if (params.genome() != null) {
      mGenomeSequences = params.genome().reader();

      try {
        mTemplateNameMap = ReaderUtils.getSequenceNameMap(mGenomeSequences);
      } catch (final IOException e) {
        throw new NoTalkbackSlimException(ErrorType.READING_ERROR, mParams.genome().toString());
      }
    } else {
      mGenomeSequences = null;
      mTemplateNameMap = null;
      mNumberTemplateSequences = 0;
      mProgressNT = 0;
      mTotalTemplateLength = 0;
    }
    mFilterParams = filteringParams;
  }

  /**
   * Perform optional one-off initialization based on the SAM header. This
   * default implementation does nothing.
   * @param header the SAM header
   * @throws java.io.IOException if there is an I/O problem
   */
  protected void init(SAMFileHeader header) throws IOException { }

  /**
   * Return true to indicate that the record was processed without errors.
   *
   * @param rec the record to process.
   * @return true if the record was processed correctly, false if the record contained errors.
   * @throws java.io.IOException if there is an I/O problem
   */
  protected abstract boolean processRecord(final SAMRecord rec) throws IOException;

  /**
   * Flush a range of template positions within the current template.
   * This is called automatically at the end of <code>exec()</code>, but
   * can also be called periodically from within <code>processRecord</code>.
   * If you do call flush, its return value should be put into
   * <code>mPreviousStart</code>, and <code>mTemplateName</code>
   * and <code>mTemplateLength</code> must be kept up to date.
   *
   * @param start first position to be output.
   * @param last one past last position to be output.
   * @return last position.
   * @exception IOException if an I/O error occurs.
   */
  public abstract int flush(final int start, final int last) throws IOException;

  /**
   * Man this interface is bad. Required because flush in variant task keeps
   * a buffer it doesn't purge until it is a certain amount past its last interesting call
   * so it needs a notify after every sequence.
   * this provides the ability to clear the final results (which isn't known in the sub class)
   * @throws java.io.IOException if there is an I/O problem
   */
  protected abstract void finalPostFlush() throws IOException;

  /**
   * Validate record
   * @param rec the record to validate
   * @return true if valid, false otherwise
   */
  protected boolean validateRecord(final SAMRecord rec) {
    final Integer nh = SamUtils.getNHOrIH(rec);
    return nh == null || nh > 0;
  }

  @Override
  protected void exec() throws IOException {
    try {
      final RecordIterator<SAMRecord> iterator;
      // mParam.thread returns T - 1 threads
      final SAMFileHeader header = SamUtils.getUberHeader(mGenomeSequences, mParams.mapped());
      if (mParams.ioThreads() < 1) {
        iterator = new MultifileIterator(new SamReadingContext(mParams.mapped(), 1, mFilterParams, header, mGenomeSequences));
      } else {
        final SingletonPopulatorFactory<SAMRecord> pf = new SingletonPopulatorFactory<>(new SamRecordPopulator());
        final SamReadingContext context = new SamReadingContext(mParams.mapped(), mParams.ioThreads(), mFilterParams, header, mGenomeSequences);
        iterator = new ThreadedMultifileIterator<>(context, pf);
      }
      try {
        if (mTemplateNameMap != null) {   // Perform header consistency check
          for (final String name : SamUtils.getSequenceNames(iterator.header())) {
            if (!mTemplateNameMap.containsKey(name)) {
              throw new NoTalkbackSlimException("SAM header sequence " + name + " is not contained in the supplied template.");
            }
          }
        }
        Diagnostic.developerLog("threads = " + mParams.ioThreads() + " : " + iterator.getClass().getName());
        init(iterator.header());
      } catch (final RuntimeException e) {
        iterator.close();
        throw e;
      }
      //setup progress
      int refStartOffset = 0;

      //System.err.println("threads = " + mParams.threads() + " : " + innerIt.getClass().getName());
      try {
        boolean first = true;
        long invalidRecords = 0;
        long validRecords = 0;
        SAMRecord rec;
        int templatesCompleted = 0;
        long templateNTCompleted = 0;
        long referenceLength = 0;
        String lastTemplateName = null;
        long nextProgressOutput = mProgressNT + refStartOffset;
        long nextRefProgressOutput = mProgressNT + refStartOffset;
        long refDivisions = mProgressNT;
        if (mGenomeSequences != null) {
          if (mFilterParams.restriction() == null) {
            mNumberTemplateSequences = mGenomeSequences.numberSequences();
          } else {
            mNumberTemplateSequences = 1;
          }
        }
        final boolean outputTotalSequences = mNumberTemplateSequences != 0;
        while (iterator.hasNext()) {
          rec = iterator.next();
          if (validateRecord(rec) && processRecord(rec)) {
            validRecords++;
          } else {
            invalidRecords++;
            if (invalidRecords <= 5) {
              Diagnostic.warning(WarningType.SAM_BAD_FORMAT_WARNING1, rec.getSAMString().trim());
              Diagnostic.userLog("Invalid record: " + rec.getSAMString().trim());
            }
          }

          //setup progress
          if (first) {
            //otherwise we can't tell what the starting point is (restriction is overlapping so start position of first alignment could be anywhere depending on the size of the read).
            if (mGenomeSequences != null) {
              if (mFilterParams.restriction() == null) {
                mTotalTemplateLength = mGenomeSequences.totalLength();
              } else {
                final SAMSequenceDictionary dict = iterator.header().getSequenceDictionary();
                refStartOffset = Math.min(mFilterParams.restriction().resolveStart() + 1, rec.getAlignmentStart());
                mTotalTemplateLength = mFilterParams.restriction().resolveEnd(dict) - refStartOffset;
              }
              mProgressNT = Math.max(mTotalTemplateLength / PROGRESS_NT_DIV, 1);
            }
            first = false;
          }

          final String refName = rec.getReferenceName();
          if (!refName.equals(lastTemplateName)) {
            if (lastTemplateName != null) {
              templatesCompleted++;
              Diagnostic.progress("Reference " + lastTemplateName + " completed (" + templatesCompleted + (outputTotalSequences ? "/" + mNumberTemplateSequences : "") + ")");
              templateNTCompleted += referenceLength;
            }
            if (mGenomeSequences != null) {
              if (mFilterParams.restriction() == null) {
                referenceLength = getSequenceLength(mGenomeSequences, mTemplateNameMap, refName);
                refDivisions = Math.max(referenceLength / PROGRESS_NT_DIV, 1);
                nextRefProgressOutput = refDivisions;
              } else {
                referenceLength = mTotalTemplateLength;
                refDivisions = Math.max(referenceLength / PROGRESS_NT_DIV, 1);
                nextRefProgressOutput = refDivisions + refStartOffset;
              }
            }
            lastTemplateName = refName;
            nextProgressOutput = mProgressNT != 0 ? mProgressNT : PROGRESS_NT;
          }
          final int alignmentStart = rec.getAlignmentStart();
          if (referenceLength != 0 && alignmentStart >= nextRefProgressOutput) {
            //To allow for gaps larger than refDivisions between records
            nextRefProgressOutput = alignmentStart - ((alignmentStart - refStartOffset) % refDivisions);
            Diagnostic.progress("Processed " + MathUtils.round((100.0 * (nextRefProgressOutput - refStartOffset)) / referenceLength) + "% of " + refName);
            nextRefProgressOutput += refDivisions;
          }
          if (alignmentStart >= nextProgressOutput) {
            final long progressNT = mProgressNT != 0 ? mProgressNT : PROGRESS_NT;
            //To allow for gaps larger than PROGRESS_NT between records
            nextProgressOutput = alignmentStart - ((alignmentStart - refStartOffset) % progressNT);
            if (mTotalTemplateLength != 0) {
              Diagnostic.progress("Processed " + MathUtils.round((100.0 * (templateNTCompleted + (nextProgressOutput - refStartOffset))) / mTotalTemplateLength) + "% of total");
            } else {
              Diagnostic.progress("Processed " + nextProgressOutput + " NT of " + refName);
            }
            nextProgressOutput += progressNT;
          }
        }
        if (mTemplateName != null) {
          flush(mPreviousStart, mTemplateLength);
          finalPostFlush();
          templatesCompleted++;
          Diagnostic.progress("Reference " + mTemplateName + " completed (" + templatesCompleted + (outputTotalSequences ? "/" + mNumberTemplateSequences : "") + ")");
        }
        iterator.close();
        invalidRecords += iterator.getInvalidRecordsCount();
        final String invalidRecordsWarning = invalidRecords + " records skipped because of SAM format problems.";
        if (invalidRecords > 0) {
          Diagnostic.warning(invalidRecordsWarning);
        } else {
          Diagnostic.userLog(invalidRecordsWarning);
        }
        if (iterator.getFilteredRecordsCount() > 0) {
          Diagnostic.userLog(iterator.getFilteredRecordsCount() + " records skipped due to input filtering criteria");
        }
        Diagnostic.userLog(validRecords + " records processed");
      } finally {
        iterator.close();
      }
    } finally {
      if (mGenomeSequences != null) {
        mGenomeSequences.close();
      }
    }
  }

  /**
   * @param genomeSequences the sequence reader
   * @param templateNameMap the template name map
   * @param name the template name
   * @return the sequences bytes
   * @throws java.io.IOException if there is an I/O problem
   */
  protected static byte[] getSequenceBytes(final SequencesReader genomeSequences, final Map<String, Long> templateNameMap, final String name) throws IOException {
    final Long sequenceId = templateNameMap.get(name);
    if (sequenceId == null) {
      return null;
    }
    return genomeSequences.read(sequenceId);
  }

  /**
   * @param genomeSequences the sequence reader
   * @param templateNameMap the template name map
   * @param name the template name
   * @return the sequences bytes
   * @throws java.io.IOException if there is an I/O problem
   */
  protected static int getSequenceLength(final SequencesReader genomeSequences, final Map<String, Long> templateNameMap, final String name) throws IOException {
    final Long sequenceId = templateNameMap.get(name);
    if (sequenceId == null) {
      return -1;
    }
    return genomeSequences.length(sequenceId);
  }
}

