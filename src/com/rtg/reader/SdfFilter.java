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
package com.rtg.reader;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.Queue;

import com.rtg.launcher.CommonFlags;
import com.rtg.mode.DNA;
import com.rtg.mode.Protein;
import com.rtg.mode.Residue;
import com.rtg.mode.SequenceType;
import com.rtg.util.Constants;
import com.rtg.util.StringUtils;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.Validator;
import com.rtg.util.diagnostic.CliDiagnosticListener;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.diagnostic.SlimException;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.LogStream;

/**
 * Performs various filtering operations on SDF files.
 *
 */
public final class SdfFilter {

  private SdfFilter() { }
  static final String APPLICATION_NAME = "rtg sdffilter";

  private static final String USAGE_HEADER1 = "-i SDF";
  private static final String USAGE_HEADER2 = "[OPTION]... -n INT -o DIR -l SDF -r SDF";
  static final String FLAG_USAGE_HEADER = USAGE_HEADER1 + StringUtils.LS + USAGE_HEADER2;

  static final String INPUT_FLAG = "input";
  static final String LEFT_INPUT_FLAG = "left-input";
  static final String RIGHT_INPUT_FLAG = "right-input";
  static final String OUTPUT_FLAG = "output";
  static final String EXPAND_FLAG = "expand-n";
  static final String ALTERNATE_FLAG = "alternate";
  static final String RC_FLAG = "reverse-complement";

  /**
   * Gets the <code>CFlags</code> object.
   *
   * @param out stdout
   * @param err stderr
   * @return the <code>CFlags</code> object
   */
  public static CFlags getCFlags(final Appendable out, final Appendable err) {
    final CFlags flags = new CFlags(APPLICATION_NAME, out, err);
    flags.registerExtendedHelp();
    flags.setRemainderHeader(FLAG_USAGE_HEADER);

    flags.registerOptional('i', INPUT_FLAG, File.class, "SDF", "input SDF directory");
    flags.registerOptional('l', LEFT_INPUT_FLAG, File.class, "SDF", "left input SDF directory");
    flags.registerOptional('r', RIGHT_INPUT_FLAG, File.class, "SDF", "right input SDF directory");
    flags.registerRequired('o', OUTPUT_FLAG, File.class, "DIR", "output base directory (must be empty if present)");

    flags.registerOptional('n', EXPAND_FLAG, Integer.class, "INT", "number of Ns to expand");

    flags.registerOptional('a', ALTERNATE_FLAG, "alternate expansion for CG data");
    flags.registerOptional('c', RC_FLAG, "reverse complement all sequences");

    flags.setValidator(VALIDATOR);
    return flags;
  }

  private static final Validator VALIDATOR = new Validator() {
    @Override
    public boolean isValid(final CFlags flags) {

      if (flags.isSet(INPUT_FLAG)) {
        if (!CommonFlags.validateSDF(flags, INPUT_FLAG)) {
          return false;
        }

        if (flags.isSet(LEFT_INPUT_FLAG) || flags.isSet(RIGHT_INPUT_FLAG)) {
          flags.setParseMessage("Cannot set --" + LEFT_INPUT_FLAG + " or --" + RIGHT_INPUT_FLAG + " at the same time as --" + INPUT_FLAG);
          return false;
        }
      } else if (flags.isSet(LEFT_INPUT_FLAG) && flags.isSet(RIGHT_INPUT_FLAG)) {
        if (!CommonFlags.validateSDF(flags, LEFT_INPUT_FLAG)) {
          return false;
        }

        if (!CommonFlags.validateSDF(flags, RIGHT_INPUT_FLAG)) {
          return false;
        }
      } else {
        flags.setParseMessage("Must set either --" + INPUT_FLAG + " or both of --" + LEFT_INPUT_FLAG + " and --" + RIGHT_INPUT_FLAG);
        return false;
      }

      if (!CommonFlags.validateOutputDirectory((File) flags.getValue(OUTPUT_FLAG))) {
        return false;
      }

      if (flags.isSet(EXPAND_FLAG)) {
        final int n = (Integer) flags.getValue(EXPAND_FLAG);
        if (n > 1) {
          flags.setParseMessage("Only 1 N can be expanded");
          return false;
        }
      }
      return true;
    }
  };

  void expand(final SequencesReader reader, final File output, final int n, final PrereadArm arm) throws IOException {
    final SequenceDataSource source;
    if (n == 1) {
      source = new Expand1SequenceDataSource(reader);
    } else {
      source = new CopyDataSource(reader);
    }
    final SequencesWriter writer = new SequencesWriter(source, output, Constants.MAX_FILE_SIZE, reader.getPrereadType(), reader.compressed());
    writer.setSdfId(reader.getSdfId());
    writer.setPrereadArm(arm);
    writer.processSequences();
    source.close();
  }

  void expandPair(final SequencesReader leftReader, final SequencesReader rightReader, final File leftOut, final File rightOut, final int n, final boolean alternate) throws IOException {
    final byte[] leftSequence = new byte[(int) leftReader.maxLength()];
    final byte[] rightSequence = new byte[(int) rightReader.maxLength()];
    final byte[] leftSequenceQuality = new byte[(int) leftReader.maxLength()];
    final byte[] rightSequenceQuality = new byte[(int) rightReader.maxLength()];

    final Residue[] values;
    if (leftReader.type().equals(SequenceType.DNA)) {
      values = DNA.values();
    } else {
      values = Protein.values();
    }

    final boolean hasQuality = leftReader.hasQualityData() && rightReader.hasQualityData();

    try (SdfWriter rightWriter = new SdfWriter(rightOut, Constants.MAX_FILE_SIZE, rightReader.getPrereadType(), hasQuality, true, rightReader.compressed(), rightReader.type());
         SdfWriter leftWriter = new SdfWriter(leftOut, Constants.MAX_FILE_SIZE, leftReader.getPrereadType(), hasQuality, true, leftReader.compressed(), leftReader.type())) {
      // Force arms so user can use filter to switch arms if desired
      leftWriter.setPrereadArm(PrereadArm.LEFT);
      leftWriter.setSdfId(leftReader.getSdfId());
      rightWriter.setPrereadArm(PrereadArm.RIGHT);
      rightWriter.setSdfId(rightReader.getSdfId());

      while (leftReader.nextSequence() && rightReader.nextSequence()) {
        final int leftLength = leftReader.readCurrent(leftSequence);
        int leftNpos = -1;
        if (n > 0) {
          for (int i = 0; i < leftLength; i++) {
            if (values[leftSequence[i]].ignore()) {
              leftNpos = i;
              break;
            }
          }
        }
        final int rightLength = rightReader.readCurrent(rightSequence);
        int rightNpos = -1;
        if (n > 0) {
          for (int i = 0; i < rightLength; i++) {
            if (values[rightSequence[i]].ignore()) {
              rightNpos = i;
              break;
            }
          }
        }
        if (hasQuality) {
          leftReader.readCurrentQuality(leftSequenceQuality);
          rightReader.readCurrentQuality(rightSequenceQuality);
        }

        if (leftNpos == -1 && rightNpos == -1) {
          writeSequence(leftWriter, leftReader.currentName(),
            leftSequence, hasQuality ? leftSequenceQuality : null, leftLength);
          writeSequence(rightWriter, rightReader.currentName(),
            rightSequence, hasQuality ? rightSequenceQuality : null, rightLength);
        } else {
          final int firstValid = leftReader.type().firstValid();

          if (alternate) {
            // preferably have a N after the first 5bp in left
            // and before last 5bp in right
            if (leftNpos != -1 && leftNpos < 5) {
              for (int i = leftNpos + 1; i < leftLength; i++) {
                if (values[leftSequence[i]].ignore()) {
                  leftNpos = i;
                  break;
                }
              }
            }
            // since we search first to last for the right N, we already have
            // found the N if it is before the last 5bp

            // expand each side (if needed) separately
            if (leftNpos != -1) {
              // remember the N that the left side had
              final byte originalN = leftSequence[leftNpos];

              for (int i = firstValid; i < values.length; i++) {
                final String suffix = "-exp-left-N->" + values[i];
                leftSequence[leftNpos] = (byte) i;

                writeSequence(leftWriter, leftReader.currentName() + suffix,
                  leftSequence, hasQuality ? leftSequenceQuality : null, leftLength);
                writeSequence(rightWriter, rightReader.currentName() + suffix,
                  rightSequence, hasQuality ? rightSequenceQuality : null, rightLength);
              }

              // restore the N to the left side
              leftSequence[leftNpos] = originalN;
            }

            if (rightNpos != -1) {
              for (int i = firstValid; i < values.length; i++) {
                final String suffix = "-exp-right-N->" + values[i];
                rightSequence[rightNpos] = (byte) i;

                writeSequence(leftWriter, leftReader.currentName() + suffix,
                  leftSequence, hasQuality ? leftSequenceQuality : null, leftLength);
                writeSequence(rightWriter, rightReader.currentName() + suffix,
                  rightSequence, hasQuality ? rightSequenceQuality : null, rightLength);
              }
            }
          } else {
            // expand both sides at once
            for (int i = firstValid; i < values.length; i++) {
              final String suffix = "-expN->" + values[i];

              if (leftNpos != -1) {
                leftSequence[leftNpos] = (byte) i;
              }
              writeSequence(leftWriter, leftReader.currentName() + suffix,
                leftSequence, hasQuality ? leftSequenceQuality : null, leftLength);

              if (rightNpos != -1) {
                rightSequence[rightNpos] = (byte) i;
              }
              writeSequence(rightWriter, rightReader.currentName() + suffix,
                rightSequence, hasQuality ? rightSequenceQuality : null, rightLength);
            }
          }
        }
      }
    }
  }

  private static void writeSequence(SdfWriter writer, String name, byte[] sequence, byte[] quality, int length) throws IOException {
    writer.startSequence(name);
    writer.write(sequence, quality, length);
    writer.endSequence();
  }

  /**
   * Runs the program and System.exits
   *
   * @param args the command line arguments
   */
  public static void main(final String[] args) {
    System.exit(mainInit(args, FileUtils.getStdoutAsOutputStream(), System.err));
  }

  /**
   * Runs the program and returns the exit code.
   *
   * @param args the command line arguments
   * @param out stdout
   * @param err stderr
   * @return 0 for success, 1 for failure
   * @throws SlimException wrapping any exceptions that occur
   */
  public static int mainInit(final String[] args, final OutputStream out, final PrintStream err) {
    return mainExec(args, out, err, PrereadVerifier.initialLog(APPLICATION_NAME));
  }

  private static SequencesReader reader(final File input, final boolean rc) throws IOException {
    try {
      return rc
        ? new ReverseComplementingReader(SequencesReaderFactory.createDefaultSequencesReaderCheckEmpty(input))
        : SequencesReaderFactory.createDefaultSequencesReaderCheckEmpty(input);
    } catch (final FileNotFoundException e) {
      throw new NoTalkbackSlimException(ErrorType.FILE_NOT_FOUND, input.toString());
    }
  }

  static int mainExec(final String[] args, final OutputStream out, final PrintStream err, final LogStream initLog) {
    CliDiagnosticListener listener = null;
    final PrintStream outStream = new PrintStream(out);
    try {
      listener = PrereadVerifier.initializeLogs(args, err, initLog);
      final CFlags flags = getCFlags(outStream, err);
      if (flags.setFlags(args)) {
        final SdfFilter filter = new SdfFilter();

        final File output = (File) flags.getValue(OUTPUT_FLAG);
        final int n;
        if (flags.isSet(EXPAND_FLAG)) {
          n = (Integer) flags.getValue(EXPAND_FLAG);
        } else {
          n = 0;
        }

        if (flags.isSet(INPUT_FLAG)) {
          try (SequencesReader input = reader((File) flags.getValue(INPUT_FLAG), flags.isSet(RC_FLAG))) {
            filter.expand(input, output, n, input.getArm());
          }
        } else {
          try (final SequencesReader leftReader = reader((File) flags.getValue(LEFT_INPUT_FLAG), flags.isSet(RC_FLAG));
               final SequencesReader rightReader = reader((File) flags.getValue(RIGHT_INPUT_FLAG), flags.isSet(RC_FLAG))) {
            FileUtils.ensureOutputDirectory(output);
            final File leftOut = new File(output, "left");
            final File rightOut = new File(output, "right");
            filter.expandPair(leftReader, rightReader, leftOut, rightOut, n, flags.isSet(ALTERNATE_FLAG));
          }
        }

        Diagnostic.deleteLog();
        return 0;
      } else {
        Diagnostic.deleteLog();
        return 1;
      }
    } catch (final SlimException e) {
      throw e;
    } catch (final Throwable t) { //catch everything except SlimException
      throw new SlimException(t);
    } finally {
      outStream.flush();
      Diagnostic.removeListener(listener);
      Diagnostic.closeLog();
    }
  }

  private static class Expand1SequenceDataSource implements SequenceDataSource {
    private final SequencesReader mReader;
    private final Residue[] mValues;
    private final byte[] mCurrentSequence;
    private final byte[] mCurrentSequenceQuality;
    private String mCurrentName;
    private final Queue<byte[]> mSequenceQueue;
    private final Queue<byte[]> mSequenceQualityQueue;
    private final Queue<String> mNameQueue;
    private long mMaxLength = Long.MIN_VALUE;
    private long mMinLength = Long.MAX_VALUE;

    public Expand1SequenceDataSource(SequencesReader reader) {
      mReader = reader;
      if (mReader.maxLength() > Integer.MAX_VALUE) {
        throw new IllegalArgumentException("Sequence length too long");
      }
      mCurrentSequence = new byte[(int) mReader.maxLength()];
      mCurrentSequenceQuality = new byte[(int) mReader.maxLength()];
      mCurrentName = null;
      if (mReader.type().equals(SequenceType.DNA)) {
        mValues = DNA.values();
      } else {
        mValues = Protein.values();
      }
      mSequenceQueue = new LinkedList<>();
      mSequenceQualityQueue = new LinkedList<>();
      mNameQueue = new LinkedList<>();
    }

    @Override
    public SequenceType type() {
      return mReader.type();
    }

    @Override
    public boolean nextSequence() throws IOException {
      if (mSequenceQueue.size() > 0) {
        byte[] b = mSequenceQueue.remove();
        System.arraycopy(b, 0, mCurrentSequence, 0, b.length);
        b = mSequenceQualityQueue.remove();
        System.arraycopy(b, 0, mCurrentSequenceQuality, 0, b.length);
        mCurrentName = mNameQueue.remove();
      } else {
        if (!mReader.nextSequence()) {
          return false;
        }
        mReader.readCurrent(mCurrentSequence);
        mReader.readCurrentQuality(mCurrentSequenceQuality);
        mCurrentName = mReader.currentName();
        for (int i = 0; i < mReader.currentLength(); i++) {
          if (mValues[mCurrentSequence[i]].ignore()) {
            mCurrentName = expand(mCurrentSequence, mCurrentSequenceQuality, mCurrentName, i, mReader.currentLength());
            break;
          }
        }
      }
      mMinLength = Math.min(mMinLength, currentLength());
      mMaxLength = Math.max(mMaxLength, currentLength());
      return true;
    }

    private String expand(final byte[] seq, final byte[] quality, final String name, final int pos, final int length) {
      final int firstValid = mReader.type().firstValid();
      final byte[] qualityCopy = Arrays.copyOf(quality, length);
      for (int i = firstValid + 1; i < mValues.length; i++) {
        final byte[] seqCopy = Arrays.copyOf(seq, length);
        seqCopy[pos] = (byte) i;
        mSequenceQueue.add(seqCopy);
        mSequenceQualityQueue.add(qualityCopy);
        mNameQueue.add(name + "-expN->" + mValues[i].toString());
      }

      seq[pos] = (byte) firstValid;
      return name + "-expN->" + mValues[firstValid].toString();
    }

    @Override
    public String name() {
      return mCurrentName;
    }

    @Override
    public boolean hasQualityData() {
      return mReader.hasQualityData();
    }

    @Override
    public void close() throws IOException {
      mReader.close();
    }

    @Override
    public void setDusting(boolean val) {
      throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public long getWarningCount() {
      return 0;
    }

    @Override
    public int currentLength() {
      return mReader.currentLength();
    }

    @Override
    public byte[] qualityData() throws IllegalStateException, IOException {
      return mCurrentSequenceQuality;
    }

    @Override
    public byte[] sequenceData() throws IllegalStateException, IOException {
      return mCurrentSequence;
    }

    @Override
    public long getDusted() {
      return 0;
    }

    @Override
    public long getMaxLength() {
      return mMaxLength;
    }

    @Override
    public long getMinLength() {
      return mMinLength;
    }
  }

  private static class CopyDataSource implements SequenceDataSource {
    private final SequencesReader mReader;
    private final byte[] mCurrentSequence;
    private final byte[] mCurrentSequenceQuality;
    private String mCurrentName;
    private long mMinLength = Long.MAX_VALUE;
    private long mMaxLength = Long.MIN_VALUE;

    public CopyDataSource(SequencesReader reader) {
      mReader = reader;
      if (mReader.maxLength() > Integer.MAX_VALUE) {
        throw new IllegalArgumentException("Sequence length too long");
      }
      mCurrentSequence = new byte[(int) mReader.maxLength()];
      mCurrentSequenceQuality = new byte[(int) mReader.maxLength()];
      mCurrentName = null;
    }

    @Override
    public SequenceType type() {
      return mReader.type();
    }

    @Override
    public boolean nextSequence() throws IOException {
      if (!mReader.nextSequence()) {
        return false;
      }
      mReader.readCurrent(mCurrentSequence);
      mReader.readCurrentQuality(mCurrentSequenceQuality);
      mCurrentName = mReader.currentName();
      mMinLength = Math.min(mMinLength, currentLength());
      mMaxLength = Math.max(mMaxLength, currentLength());
      return true;
    }

    @Override
    public String name() {
      return mCurrentName;
    }

    @Override
    public byte[] sequenceData() {
      return mCurrentSequence;
    }

    @Override
    public byte[] qualityData() {
      return mCurrentSequenceQuality;
    }

    @Override
    public boolean hasQualityData() {
      return mReader.hasQualityData();
    }

    @Override
    public void close() throws IOException {
      mReader.close();
    }

    @Override
    public void setDusting(boolean val) {
      throw new UnsupportedOperationException();
    }

    @Override
    public int currentLength() {
      return mReader.currentLength();
    }

    @Override
    public long getWarningCount() {
      return 0;
    }

    @Override
    public long getDusted() {
      return 0;
    }

    @Override
    public long getMaxLength() {
      return mMaxLength;
    }

    @Override
    public long getMinLength() {
      return mMinLength;
    }
  }

}
