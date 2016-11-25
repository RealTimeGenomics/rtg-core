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
package com.rtg.variant.sv;

import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStream;

import com.rtg.launcher.CommonFlags;
import com.rtg.tabix.TabixIndexer;
import com.rtg.tabix.UnindexableDataException;
import com.rtg.util.StringUtils;
import com.rtg.util.Utils;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;


/**
 */
public class SvInterestingRegionExtractor implements Closeable {

  private static final byte[] TAB_BYTES = StringUtils.TAB.getBytes();

  private static final String FILENAME_INTERESTING = "sv_interesting.bed";

  private int mCurrentMaxValuePos = 0;
  private int mCurrentRegionStart = 0;
  private double mCurrentRegionMaxScore = 0;
  private double mCurrentRegionTotalScore = 0;
  private int mCurrentRegionPositions = 0;
  private int mRegionChanges = 0;

  private int mCurrentTemplateLength = 0;

  private static final int NORMAL_POS = 0;  //at this stage, assume NORMAL is at index 0

  private final InterestingRegionEmitter mEmitter;

  SvInterestingRegionExtractor(SvToolParams params) throws IOException {
    mEmitter = new InterestingRegionEmitter(params.outStream(FILENAME_INTERESTING), 0);
  }

  SvInterestingRegionExtractor(OutputStream outputStream, Integer normalRegionLength) {
    mEmitter = new InterestingRegionEmitter(outputStream, normalRegionLength);
  }

  /**
   * Call this each time the template changes, and once at the end to clean up last region.
   * @param templateName name of the new template
   * @param templateLength length of the new template
   * @throws IOException if an IO exception occurs when writing to the stream
   */
  void setTemplate(String templateName, int templateLength) throws IOException {
    if (mEmitter.mTemplateName != null) {
      if (mCurrentRegionPositions > 0) {
        mEmitter.mergeOrOutputRegion(mCurrentRegionStart, mCurrentTemplateLength, mRegionChanges, mCurrentRegionMaxScore, mCurrentRegionTotalScore, mCurrentRegionPositions);
        mCurrentRegionMaxScore = 0;
        mCurrentRegionPositions = 0;
        mCurrentRegionTotalScore = 0;
        mRegionChanges = 0;
        mCurrentMaxValuePos = NORMAL_POS;
      }
      mEmitter.outputRegion();
    }
    mEmitter.mTemplateName = templateName;
    mCurrentTemplateLength = templateLength;
  }

  private static class InterestingRegionEmitter implements Closeable {
    private int mStartPos;
    private int mEndPos = -1;
    private int mChanges;
    private double mMaxScore;
    private double mTotalScore;
    private int mRegions;

    private final int mMinNormalRegionLength;

    String mTemplateName = null;

    private boolean mHeaderEmitted = false;

    private OutputStream mOutputStream = null;


    InterestingRegionEmitter(OutputStream outputStream, Integer normalRegionLength) {
      mOutputStream = outputStream;
      if (normalRegionLength != null) {
        mMinNormalRegionLength = normalRegionLength;
      } else {
        mMinNormalRegionLength = 0;
      }
    }

    void mergeOrOutputRegion(int oneBasedStartPos, int oneBasedEndPos, int changes, double maxScore, double totalScore, int regions) throws IOException {
      if (mEndPos != -1 && oneBasedStartPos - mEndPos < mMinNormalRegionLength) {
        //merge the regions
        mEndPos = oneBasedEndPos;
        mChanges += changes + 2;        //additional 1 for the normality region in between
        mMaxScore = Math.max(mMaxScore, maxScore);
        mTotalScore += totalScore;
        mRegions += regions;
      } else {
        //output the first region, set values to new ones.
        outputRegion();
        mStartPos = oneBasedStartPos;
        mEndPos = oneBasedEndPos;
        mChanges = changes;
        mMaxScore = maxScore;
        mTotalScore = totalScore;
        mRegions = regions;

      }
    }


    void outputHeader() throws IOException {
      mOutputStream.write("#chr\t".getBytes());
      mOutputStream.write("start\t".getBytes());
      mOutputStream.write("end\t".getBytes());
      mOutputStream.write("areas\t".getBytes());
      mOutputStream.write("maxscore\t".getBytes());
      mOutputStream.write("average".getBytes());
      mOutputStream.write(StringUtils.LS.getBytes());
    }


    void outputRegion() throws IOException {
      if (mEndPos == -1) {      //if no current region info, just ignore.
        return;
      }
      if (!mHeaderEmitted) {
        outputHeader();
        mHeaderEmitted = true;
      }

      //Bed is 0-based positions, so convert start/end positions
      mOutputStream.write(mTemplateName.getBytes());
      mOutputStream.write(TAB_BYTES);
      mOutputStream.write(("" + (mStartPos - 1)).getBytes());
      mOutputStream.write(TAB_BYTES);
      mOutputStream.write(("" + Math.max(0, mEndPos - 1)).getBytes());
      mOutputStream.write(TAB_BYTES);
      mOutputStream.write(("" + (mChanges + 1)).getBytes());
      mOutputStream.write(TAB_BYTES);
      mOutputStream.write(Utils.realFormat(mMaxScore, 4).getBytes());
      mOutputStream.write(TAB_BYTES);
      mOutputStream.write(Utils.realFormat(mTotalScore / mRegions, 4).getBytes());
      mOutputStream.write(StringUtils.LS.getBytes());
      clear();
    }

    private void clear() {
      mStartPos = 0;
      mEndPos = -1;
      mChanges = 0;
      mMaxScore = 0;
      mTotalScore = 0;
      mRegions = 0;
    }

    @Override
    public void close() throws IOException {
      mOutputStream.close();
    }

  }


  void processValue(int oneBasedPosition, double[] values, int maxValuePos) throws IOException {
    if (maxValuePos != NORMAL_POS) {    //this position is in an interesting region
      if (mCurrentMaxValuePos == NORMAL_POS) {  // in fact, this is the start of an interesting region
        mCurrentRegionStart = oneBasedPosition;
        mCurrentRegionMaxScore = 0;
        mCurrentRegionPositions = 0;
        mCurrentRegionTotalScore = 0;
        mRegionChanges = 0;
      } else if (mCurrentMaxValuePos != maxValuePos) {
        mRegionChanges++;
      }
      mCurrentMaxValuePos = maxValuePos;
      final double thisMaxValue = values[maxValuePos];
      if (thisMaxValue > mCurrentRegionMaxScore) {
        mCurrentRegionMaxScore = thisMaxValue;
      }
      mCurrentRegionPositions++;
      mCurrentRegionTotalScore += thisMaxValue;

    } else {
      if (mCurrentMaxValuePos != NORMAL_POS) { //but it was the end of an interesting region
        //outputRegion(position);
        mEmitter.mergeOrOutputRegion(mCurrentRegionStart, oneBasedPosition, mRegionChanges, mCurrentRegionMaxScore, mCurrentRegionTotalScore, mCurrentRegionPositions);
        mCurrentRegionPositions = 0;
      }
      mCurrentMaxValuePos = NORMAL_POS;
    }
  }

  @Override
  public void close() throws IOException {
    mEmitter.close();
  }


  private static void initFlags(CFlags flags) {
    flags.registerRequired('i', "input", File.class, "FILE", "SV file to process").setCategory(INPUT_OUTPUT);
    flags.registerOptional('o', CommonFlags.OUTPUT_FLAG, File.class, "FILE", "output BED file").setCategory(INPUT_OUTPUT);

    //set length of normal region which must be present after a breakpoint (for merging close regions)
    flags.registerOptional('m', "merge-regions", Integer.class, "INT", "number of windows which regions must be apart to be considered separate").setCategory(CommonFlagCategories.FILTERING);
  }

  /**
   * Extracts interesting regions from a pre-existing Sv bayesian output file.
   * @param args the arguments
   * @throws IOException if an IO exception occurs
   */
  public static void main(String[] args) throws IOException {
    final CFlags flags = new CFlags();
    initFlags(flags);
    if (!flags.setFlags(args)) {
      return;
    }

    final File input = (File) flags.getValue("input");
    final File output = (File) flags.getValue(CommonFlags.OUTPUT_FLAG);
    final Integer merge = (Integer) flags.getValue("merge-regions");

    final boolean gzipOutput = output != null && FileUtils.isGzipFilename(output);

    int skipped = 0;
    try (SvInterestingRegionExtractor re = new SvInterestingRegionExtractor(output != null ? FileUtils.createOutputStream(output, gzipOutput, false) : FileUtils.getStdoutAsOutputStream(), merge)) {
      try (BufferedReader br = new BufferedReader(new InputStreamReader(FileUtils.createInputStream(input, false)))) {
        String line;
        double[] values = null;
        while ((line = br.readLine()) != null) {
          if (line.startsWith("#") || line.trim().length() == 0) {
            continue;
          }
          final String[] parts = line.split("\t");
          if (parts.length < 4) {
            skipped++;
            continue;
          }
          if (values == null) {
            values = new double[parts.length - 3];      // chr, pos, maxValueIndex are stripped
          }
          final String chr = parts[0];
          //if chr has changed, emit last region
          if (!chr.equals(re.mEmitter.mTemplateName)) {
            re.setTemplate(chr, 0);
          }

          final Integer pos = Integer.parseInt(parts[1]);
          final Integer maxValueIndex = Integer.parseInt(parts[parts.length - 1]);
          for (int i = 2; i < parts.length - 2; i++) {
            values[i - 2] = Double.parseDouble(parts[i]);
          }
          re.processValue(pos, values, maxValueIndex);
        }
        //need to emit the last interesting region too, if there is one
        re.setTemplate(null, -1);
      }
    }
    if (skipped > 0) {
      System.out.println(skipped + " lines were skipped for being too short.   ");
    }
    if (output != null && gzipOutput) {
      try {
        new TabixIndexer(output, TabixIndexer.indexFileName(output)).saveBedIndex();
      } catch (UnindexableDataException e) {
        Diagnostic.warning("Cannot produce TABIX index for: " + output + ": " + e.getMessage());
      }
    }
  }
}
