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
package com.rtg.index.hash.ngs;

import java.io.IOException;

import com.rtg.launcher.HashingRegion;
import com.rtg.launcher.ISequenceParams;
import com.rtg.mode.SequenceMode;
import com.rtg.mode.UnidirectionalFrame;
import com.rtg.reader.SequencesReader;
import com.rtg.util.IORunnable;
import com.rtg.util.ProgramState;
import com.rtg.util.SimpleThreadPool;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.WarningType;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 */
public class NgsHashLoopImpl extends IntegralAbstract implements NgsHashLoop {

  private static final boolean DAVE_N_HACK = false;

  private static final int WRONG_LENGTH_REPORT_LIMIT = 10;

  private static final int MASKED_CALL_PADDING_ADJUSTMENT = 64;
  /** Timer for delays in doing I/O. */
  //  public static final Timer READ_DELAY = new Timer("Read_delay");
  private final boolean mProgress;
  private final long mNumberReads;
  private boolean mReadSequencesDefined = false;
  private final long mReadProgressMask;
  private final long mTemplateProgressMask;
  private long mThreadPadding;
  /** Value to use for the generation of a valid base when an N is seen. */
  private int mUnknownVictim = 0;

  protected long mMinChunkSize = HashingRegion.DEFAULT_MIN_CHUNK_SIZE;

  /**
   * set the number of bases to pad each threads region to handle matches near the boundaries
   * @param threadPadding the number of bases to pad the thread
   */
  public void setThreadPadding(final long threadPadding) {
    this.mThreadPadding = threadPadding;
  }

  /**
   * @param numberReads total number of reads used for building.
   * @param progress true iff progress messages are to be generated.
   */
  public NgsHashLoopImpl(final long numberReads, final boolean progress) {
    this(numberReads, progress, 0x3FFFFL, 0x1FFFFL);
  }

  /**
   * @param numberReads total number of reads used for building.
   * @param progress true iff progress messages are to be generated.
   * @param readProgressMask sets the frequency of read progress (intended for testing).
   * @param templateProgressMask sets the frequency of template progress (intended for testing).
   */
  public NgsHashLoopImpl(final long numberReads, final boolean progress, final long readProgressMask, final long templateProgressMask) {
    mNumberReads = numberReads;
    mProgress = progress;
    mReadProgressMask = readProgressMask;
    mTemplateProgressMask = templateProgressMask;
  }

  @Override
  public long readLoop(final ISequenceParams params, final ReadHashFunction hashFunction, final ReadEncoder encoder, final boolean reverse) throws IOException {
    final SequenceMode mode = SequenceMode.UNIDIRECTIONAL;
    final UnidirectionalFrame frame = UnidirectionalFrame.FORWARD;

    final SequencesReader reader = params.reader();
    final long start = 0;
    final long end = reader.numberSequences();
    if (end > Integer.MAX_VALUE) {
      throw new RuntimeException("Too many reads");
    }

    assert mode.codeType().firstValid() == 1;
    final byte[] byteBuffer = makeBuffer(reader);
    int badLengthCount = 0;
    long totalLength = 0;
    for (int seq = (int) start; seq < end; ++seq) {
      mUnknownVictim = -1;
      final int readId = encoder.encode(seq);
      final int id2 = encoder.encode((int) (seq - start));
      //System.err.println("id2=" + id2);
      if ((id2 & mReadProgressMask) == 0) {
        ProgramState.checkAbort();
      }
      //System.err.println("seq=" + seq);
      final int currentLength = reader.length(seq);
      if (currentLength != hashFunction.readLength()) {
        if (badLengthCount++ < WRONG_LENGTH_REPORT_LIMIT) {
          Diagnostic.warning(WarningType.INCORRECT_LENGTH, reader.names() != null ? reader.name(seq) : ("" + seq), currentLength + "", hashFunction.readLength() + "");
        }
        hashFunction.setValues(id2, false);
        continue;
      }
      final int length = reader.read(seq, byteBuffer);
      totalLength += length;
      //System.err.println(Arrays.toString(byteBuffer));
      hashFunction.reset();
      int prev = 0;
      for (int j = 0; j < currentLength; ++j) {
        //System.err.println("j=" + j);
        final byte b = frame.code(byteBuffer, length, j);
        final int c = b - 1;
        //System.err.println("b=" + b + " c=" + c);
        //System.err.println("c=" + c);
        final byte v;
        if (c < 0) {
          if (DAVE_N_HACK || mUnknownVictim == -1) {
            mUnknownVictim = prev;
          }
          v = (byte) mUnknownVictim;
          if (!DAVE_N_HACK) {
            ++mUnknownVictim;
            mUnknownVictim &= 3;
          }
        } else {
          v = (byte) c;
          prev = c;
        }
        hashFunction.hashStep(v);
      }
      hashFunction.readAll(readId, reverse);
      hashFunction.setValues(id2, reverse);
      hashFunction.reset();
    }
    if (badLengthCount >= WRONG_LENGTH_REPORT_LIMIT) {
      Diagnostic.warning(WarningType.NUMBER_OF_INCORRECT_LENGTH, String.valueOf(badLengthCount));
    }

    mReadSequencesDefined = true;
    assert Exam.globalIntegrity(hashFunction);
    return totalLength;
  }

  /**
   * Scans an array of hashes and makes calls to
   * <code>hashFunction</code>.
   * (used when doing matched exclusion in PCR).
   * @param hashes array of hash values to be re-encoded for use by matching.
   * @param hashFunction it is given the succesive codes and does the work from there on.
   * @throws IOException If an I/O error occurs
   */
  public void pcrLoop(final long[] hashes, final ReadHashFunction hashFunction) throws IOException {
    hashFunction.setReadSequences(mNumberReads);

    final long numberSeq = hashes.length;
    final int readLength = hashFunction.readLength();
    final int shift0 = 2 * (readLength - 1);
    for (int seq = 0, id2 = 0; seq < numberSeq; ++seq, id2 += 1) {
      if ((id2 & mReadProgressMask) == 0) {
        ProgramState.checkAbort();
      }
      hashFunction.reset();
      final long hash = hashes[seq];
      //System.err.println("seq=" + seq + " hash=" + Utils.toBits(hash));
      for (int s = shift0; s >= 0; s -= 2) {
        //System.err.println("j=" + j);
        final byte c = (byte) ((hash >> s) & 3);
        //System.err.println("c=" + c);
        hashFunction.hashStep(c);
      }
      hashFunction.readAll(seq, false);
      hashFunction.setValues(id2, false);
      hashFunction.reset();
    }
    mReadSequencesDefined = true;
  }

  /**
   * Scans the buffer and makes calls to
   * <code>hashStepTemplate</code> in the hash function.
   * @param params specifies start, end and reader.
   * @param hashFunction does the work as each code is added from the template.
   * @throws IOException if an I/O error occurs.
   */
  @Override
  public void templateLoop(final ISequenceParams params, final TemplateHashFunction hashFunction) throws IOException {
    if (!mReadSequencesDefined) {
      throw new RuntimeException("Read sequences not defined");
    }
    final HashingRegion region = params.region();
    final SequencesReader reader = params.reader();
    final long start;
    final long end;
    if (region != HashingRegion.NONE) {
      start = region.getStart();
      end = region.getEnd();
    } else {
      start = 0;
      end = reader.numberSequences();
    }

    final SequenceMode mode = SequenceMode.BIDIRECTIONAL;
    final byte[] byteBuffer = makeBuffer(reader);
    //number of codes to be ignored (unknown etc).
    assert mode.codeType().firstValid() == 1;
    //the various possible frames (direction and phase for translation)
    assert mode.allFrames().length == 2;
    long nt = 0;
    long ntbatch = 0;
    long batchstart = System.currentTimeMillis();
    for (long templateId = start; templateId < end; ++templateId) {
      final int length = reader.read(templateId, byteBuffer);
      hashFunction.reset();
      hashFunction.templateSet(templateId, length);
      for (int endPosition = 0; endPosition < length; ++endPosition, ++nt) {
        if ((nt & mTemplateProgressMask) == 0) {
          ProgramState.checkAbort();
        }
        if (mProgress) {
          // every 10,000,000nt print out a debug line
          if (ntbatch >= 10000000) {
            final long bs = System.currentTimeMillis();
            final double secs = (bs - batchstart) / 1000.0;
            final long ntpersec = (int) (ntbatch / secs);
            // estimate time for human in hours
            final long threebillestimate = (int) (3000000000L / ntpersec) / 3600;
            Diagnostic.userLog("Throughput " + ntbatch + "nt in " + secs + "s, rate=" + ntpersec + " nt/s: estimate for 3 billion nt=" + threebillestimate + " hours");
            ntbatch = 0;
            batchstart = bs;
          }
          ++ntbatch;
        }
        //final byte b = frame.code(byteBuffer, length, endPosition);
        final byte b = byteBuffer[endPosition];
        final int c = b - 1;
        if (c < 0) {
          hashFunction.reset();
          continue;
        }
        hashFunction.hashStep((byte) c);
        hashFunction.templateBidirectional(endPosition);
      } //sequence
      for (int endPosition = length; endPosition < length + 64 - hashFunction.readLength(); ++endPosition) {
        hashFunction.hashStep((byte) 0);
        hashFunction.templateReverse(endPosition);
      }
      hashFunction.endSequence();
    }
    hashFunction.logStatistics();
  }
  /** Used to compute the maximum number of sequences before multiple sequences are allocated to each thread. */
  public static final int MAX_SEQUENCES = 10;


  /**
   * Scans the buffer and makes calls to
   * <code>hashStepTemplate</code> in the hash function.
   * Use threads to run this in parallel
   * @param params specifies start, end and reader.
   * @param hf model used to create on hash function for each thread.
   * @param numberThreads the number of threads to use.
   * @param threadMultiplier the number of chunks of work to allocate to each thread.
   * @throws IOException if an I/O error occurs.
   */
  @Override
  public void templateLoopMultiCore(final ISequenceParams params, final NgsHashFunction hf, final int numberThreads, int threadMultiplier) throws IOException {
    if (!mReadSequencesDefined) {
      throw new RuntimeException("Read sequences not defined");
    }
    final HashingRegion region = params.region();
    final long start;
    final long end;
    if (region != HashingRegion.NONE) {
      start = region.getStart();
      end = region.getEnd();
    } else {
      start = 0;
      end = params.reader().numberSequences();
    }
    final SequenceMode mode = SequenceMode.BIDIRECTIONAL;
    //number of codes to be ignored (unknown etc)
    assert mode.codeType().firstValid() == 1;
    assert mode.allFrames().length == 2;
    final SequencesReader reader0 = params.reader();
    final long t0 = System.currentTimeMillis();
    final SimpleThreadPool pool = new SimpleThreadPool(numberThreads, "Search", true);

    // For a 50m Yoruba run taking 40 minutes for 16 threads this is 640 cpu minutes.
    // For Threads=16 * 2 = 32 chunks, this is 20 minutes for each chunk. Trying to reduce the maximum
    // latency to Threads * 8 = could result in a maximum delay of 5 minutes for each chunk.
    // In the case of a single template sequence, the region splitting has a minimum chunk size so
    // this will avoid creating too many chunks.
    final HashingRegion[] ranges = HashingRegion.splitWorkload(reader0, params.sex(), start, end, numberThreads * threadMultiplier, mMinChunkSize, mThreadPadding);
    pool.enableBasicProgress(ranges.length);
    for (int i = 0; i < ranges.length; ++i) {
      schedule(params, hf, t0, pool, i, ranges[i]);
    }
    timeLog(t0, "parent", "Terminating", HashingRegion.NONE);
    pool.terminate();
    timeLog(t0, "parent", "Finished", HashingRegion.NONE);
  }

  private void schedule(final ISequenceParams params, final NgsHashFunction hf, final long t0, final SimpleThreadPool pool, final int i, final HashingRegion region) {
    final String name = Integer.toString(i);
    timeLog(t0, name, "Scheduling", region);
    pool.execute(new SequenceLoop(this, params, hf, region, name, t0));
  }

  private static void timeLog(final long t0, final String name, final String label, final HashingRegion region) {
    Diagnostic.userLog("Thread Search " + name + " " + label + " " + region + " " + (System.currentTimeMillis() - t0) / 1000L + " s");
  }

  static class SequenceLoop implements IORunnable {

    final NgsHashLoopImpl mParent;
    final ISequenceParams mParams;
    final NgsHashFunction mParentFunc;
    final HashingRegion mRegion;
    final String mName;
    long mT0;

    SequenceLoop(final NgsHashLoopImpl parent, final ISequenceParams params, final NgsHashFunction hf, final HashingRegion region, final String name, final long t0) {
      mParent = parent;
      mParams = params;
      mParentFunc = hf;
      mRegion = region;
      mName = name;
      mT0 = t0;
    }

    private long getPadding() {
      return mParent.mThreadPadding;
    }

    @Override
    public void run() throws IOException {
      final NgsHashFunction func = mParentFunc.threadClone(mRegion);
      timeLog(mT0, mName, "Start", mRegion);
      mT0 = System.currentTimeMillis();            // Now reset this to represent thread execution start time.
      final SequencesReader reader = mParams.subSequence(mRegion).reader();
      final long padding = getPadding();
      final byte[] byteBuffer = makeBuffer(reader, mRegion, padding);
      for (long templateId = mRegion.getStart(); templateId <= mRegion.getEnd(); ++templateId) {
        final int length = reader.length(templateId); //reader.readCurrent(byteBuffer);
        final int startPos = mRegion.getReferenceStart(templateId, padding);
        final int endPos = mRegion.getReferenceEnd(templateId, padding, length);
        final int adjStart;
        if (startPos - MASKED_CALL_PADDING_ADJUSTMENT < 0) {
          adjStart = 0;
        } else {
          adjStart = startPos - MASKED_CALL_PADDING_ADJUSTMENT;
        }
        final int adjEnd;
        if (endPos + MASKED_CALL_PADDING_ADJUSTMENT > length) {
          adjEnd = length;
        } else {
          adjEnd = endPos + MASKED_CALL_PADDING_ADJUSTMENT;
        }
        //final int measuredLength = reader.readCurrent(byteBuffer, startPos, endPos - startPos);
        final int measuredLength = reader.read(templateId, byteBuffer, adjStart, adjEnd - adjStart);

        //mParent.multiCoreFrameLoop(templateId, func, byteBuffer, measuredLength, startPos, endPos);
        mParent.multiCoreFrameLoop(templateId, func, byteBuffer, measuredLength, adjStart, adjEnd, startPos, endPos);
      }
      reader.close();
      func.threadFinish();
      timeLog(mT0, mName, "Finish", mRegion);
      func.logStatistics();
    }
  }


  private static boolean updateHashStep(TemplateHashFunction function, byte[] byteBuffer, int index) {
    final byte b = byteBuffer[index];
    final int c = b - 1;
    if (c < 0) {
      function.reset();
      return false;
    }

    function.hashStep((byte) c);
    return true;
  }

  //This was marked as broken for cg. All signs currently point to cg working in multicore -KG
  private void multiCoreFrameLoop(final long templateId, final TemplateHashFunction hashFunction, final byte[] byteBuffer, final int length, final int start, final int end, int callStart, int callEnd) throws IOException {
    //as near as I can tell "length" here is never actually used by any current implementation, so the ambiguity about whether it should be length of padded chunk or
    //length of current template sequence is unknown and currently unimportant. We are currently passing the size of the byteBuffer in use.
    hashFunction.reset();
    hashFunction.templateSet(templateId, length);
    long nt = 0;
    /* if (start > 0) {
      for (int i = (int) Math.max(start - hashFunction.readLength(), 0L); i < start; ++i) {
        final byte b = byteBuffer[i];
        final int c = b - 1;
        if (c < 0) {
          hashFunction.reset();
          continue;
        }
        hashFunction.hashStep((byte) c);
      }
    }*/
    //preload hash bits
    for (int endPosition = start; endPosition < callStart; ++endPosition, ++nt) {
      if ((nt & mTemplateProgressMask) == 0) {
        ProgramState.checkAbort();
      }
      updateHashStep(hashFunction, byteBuffer, endPosition - start);
    }
    for (int endPosition = callStart; endPosition < callEnd; ++endPosition, ++nt) {
      if ((nt & mTemplateProgressMask) == 0) {
        ProgramState.checkAbort();
      }
      //final byte b = frame.code(byteBuffer, length, endPosition);
      if (updateHashStep(hashFunction, byteBuffer, endPosition - start)) {
        hashFunction.templateBidirectional(endPosition);
      }
    }
    // flush extra positions at the end
    for (int i = 0; i < 64 - hashFunction.readLength(); ++i) {
      if (callEnd + i < end) {
        updateHashStep(hashFunction, byteBuffer, callEnd + i - start);
      } else {
        hashFunction.hashStep();
      }
      hashFunction.templateReverse(callEnd + i);
    }
    hashFunction.endSequence();
  }

  static byte[] makeBuffer(final SequencesReader reader) throws IOException {
    return makeBuffer(reader, HashingRegion.NONE, 0);
  }

  static byte[] makeBuffer(SequencesReader reader, HashingRegion r, long padding) throws IOException {
    final int longestSubSeq = (int) (r.longestSubSequence(reader) + 2 * padding);
    final int size;
    if ((long) longestSubSeq + 2 * MASKED_CALL_PADDING_ADJUSTMENT > Integer.MAX_VALUE) {
      size = Integer.MAX_VALUE;
    } else {
      size = longestSubSeq + 2 * MASKED_CALL_PADDING_ADJUSTMENT;
    }
    return new byte[size];
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mNumberReads >= 0);
    Exam.assertTrue(((mReadProgressMask + 1) & mReadProgressMask) == 0L);
    Exam.assertTrue(((mTemplateProgressMask + 1) & mTemplateProgressMask) == 0L);
    return true;
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("HashLoop");
  }


}
