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
import java.io.InputStream;
import java.util.Iterator;
import java.util.List;

import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.InformationType;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.FileUtils;

/**
 * Iterates through <code>FileStreams</code>
 */
class FileStreamIterator implements Iterator<InputStream> {

  private final Iterator<File> mIt;
  private final String mProcessLabel;
  private int mCounter = 0;
  private final int mMaxCount;
  private InputStream mLast;
  private InputStream mNext;

  private File mNextFile;
  private File mLastFile;

  /**
   * Constructs the iterator
   */
  FileStreamIterator(final List<File> files, PrereadArm arm) {
    mIt = files.iterator();
    mMaxCount = files.size();
    mLast = null;
    mNext = null;
    if (arm == null || arm == PrereadArm.UNKNOWN) {
      mProcessLabel = "";
    } else {
      mProcessLabel = (arm == PrereadArm.LEFT ? "left" : "right") + " arm ";
    }
  }

  /**
   */
  @Override
  public boolean hasNext() {
    if (mNext != null) {
      return true;
    }
    if (mIt.hasNext()) {
      mNextFile = mIt.next();
      mCounter++;
      try {
        mNext = FileUtils.createInputStream(mNextFile, true);
        if (mMaxCount == 1) {
          Diagnostic.info("Processing " + mProcessLabel + "\"" + mNextFile.toString() + "\"");
        } else {
          Diagnostic.info(InformationType.PROCESSING_ITEM_N_OF_N, true, mProcessLabel, mNextFile.toString(), Integer.toString(mCounter), Integer.toString(mMaxCount));
        }
      } catch (final FileNotFoundException fnfe) {
        throw new NoTalkbackSlimException("The file: \"" + mNextFile.getPath() + "\" either could not be found or could not be opened. The underlying error message is: \"" + fnfe.getMessage() + "\"");
      } catch (final IOException ex) {
        throw new NoTalkbackSlimException(ErrorType.IO_ERROR, "The file: \"" + mNextFile.getPath() + "\" had a problem while reading. The underlying error message is: \"" + ex.getMessage() + "\"");
      }
      return true;
    }
    return false;
  }

  /**
   */
  @Override
  public InputStream next() {
    if (mLast != null) {
      try {
        mLast.close();
      } catch (final IOException e) {
        //not the end of world, log it
        Diagnostic.userLog("Failed to close inputstream");
      }
    }
    if (mNext == null) {
      hasNext();
    }
    mLast = mNext;
    mLastFile = mNextFile;
    mNext = null;
    return mLast;
  }

  /**
   * Gets the file returned be the last call to <code>next</code>
   * @return current file, or null.
   */
  public File currentFile() {
    return mLastFile;
  }

  /**
   * Not implemented
   */
  @Override
  public void remove() {
    throw new UnsupportedOperationException();
  }
}
