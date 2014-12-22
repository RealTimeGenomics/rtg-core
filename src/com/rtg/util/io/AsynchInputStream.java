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
package com.rtg.util.io;

import java.io.IOException;
import java.io.InputStream;

/**
 * An input stream that does the reading in a separate thread.  It
 * also does buffering, during the input and also within the pipe used
 * to connect the thread to this thread.
 *
 * WARNING: Because the helper thread starts in the constructor,
 * subclasses can not contain any state themselves (as these will not
 * have been initialized before the thread starts)
 *
 */
public class AsynchInputStream extends InputStream {

  /**
   * The size of the buffer/pipe between the decompression thread and the
   * main thread.  Benchmarking showed it is best to make this fairly
   * large, probably because it allows the decompression thread to get ahead.
   */
  public static final int DEFAULT_PIPE_SIZE = 1024 * 1024;

  /** Modern disks have about 64 Kb blocks, so make it that size? */
  public static final int DEFAULT_INPUT_BUFFER_SIZE = 65536;

  /** This is package-level protection, just for testing purposes. */
  final ConcurrentByteQueue mQueue;

  private final AsynchInput mAsynchInput;

  private boolean mSeenEof = false;

  /** This is package-level protection, just for testing purposes. */
  final Thread mThread;

  /**
   * Create an asynchronous input stream to read the given file.
   *
   * @param input the input stream
   * @param pipeSize the size of the buffer between the threads.  At least 1 Kb.
   * @param bufferSize the buffer size of the input reading object.
   */
  public AsynchInputStream(InputStream input, int pipeSize, int bufferSize) {
    //Diagnostic.developerLog("new AsynchInputStream(" + pipeSize + ", " + bufferSize + ")");
    assert pipeSize >= 1024;
    mQueue = new ConcurrentByteQueue(pipeSize);
    if (input == null) {
      throw new IllegalArgumentException("File cannot be null");
    }
    mAsynchInput = new AsynchInput(input, mQueue, bufferSize);
    mThread = new Thread(mAsynchInput, "AsynchInputStream");
    mThread.setDaemon(true);
    mThread.start();
  }

  /**
   * Create an asynchronous input stream to read the given file.
   *
   * @param input the input stream
   */
  public AsynchInputStream(InputStream input) {
    this(input, DEFAULT_PIPE_SIZE, DEFAULT_INPUT_BUFFER_SIZE);
  }

  /**
   * @see java.io.InputStream#markSupported()
   * @return false
   */
  @Override
  public boolean markSupported() {
    return false;
  }

  /**
   * This also IO errors that happened in the  thread.
   */
  @Override
  public void close() throws IOException {
    if (!mSeenEof) {
      // Diagnostic.developerLog("AsynchInputStream close called before reading finished.");
      mThread.interrupt();
      try {
        mThread.join();
      } catch (final InterruptedException e) {
        throw new IOException("AsynchInputStream interrupted during close");
      } finally {
        super.close();
      }
    }
    checkException();
  }

  @Override
  public int available() {
    return mQueue.available();
  }

  private void checkException() throws IOException {
    if (mAsynchInput.mException != null) {
      try {
        throw mAsynchInput.mException;
      } finally {
        mAsynchInput.mException = null;
      }
    }
  }

  /**
   */
  @Override
  public int read(byte[] buf, int off, int len) throws IOException {
    checkException();
    try {
      final int result = mQueue.read(buf, off, len);
      if (result < 0) {
        //A result of -1 can also mean that the underlying input stream is sitting on an exception
        checkException();
        mSeenEof = true;
      }
      return result;
    } catch (final InterruptedException e) {
      throw new IOException("AsynchInputStream interrupted during read/3");
    }
  }

  /**
   */
  @Override
  public int read() throws IOException {
    throw new UnsupportedOperationException();
  }


  /**
   * This class does the reading and decompression of the input file.
   */
  private static class AsynchInput implements Runnable {

    private final InputStream mInput;

    /** the queue/pipe used to send bytes to the parent process */
    private final ConcurrentByteQueue mQueue;

    /** A buffer for the decompressed input. */
    private final byte[] mBuffer;

    volatile IOException mException = null; // tell the parent about an error.

    AsynchInput(InputStream input, ConcurrentByteQueue queue, int inputSize) {
      mInput = input;
      mQueue = queue;
      mBuffer = new byte[inputSize];
    }

    @Override
    public void run() {
      final Thread self = Thread.currentThread();
      try {
        while (!self.isInterrupted()) {
          final int size = mInput.read(mBuffer);
          if (size > 0) {
            mQueue.write(mBuffer, 0, size);
          }
          if (size <= 0) {
            break;
          }
        }
      } catch (final IOException e) {
        mException = e; // tell the other end of the pipe about this error.
      } catch (final InterruptedException e) {
        // okay, because we assume that the other thread was asking us to stop early
      } finally {
        try {
          mInput.close();
        } catch (final IOException e2) {
          if (mException == null) {
            mException = e2;
          }
        } finally {
          mQueue.close();
        }
      }
    }
  }
}

