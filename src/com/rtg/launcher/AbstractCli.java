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
package com.rtg.launcher;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;

import com.rtg.usage.UsageMetric;
import com.rtg.usage.UsageTracking;
import com.rtg.util.Constants;
import com.rtg.util.License;
import com.rtg.util.StringUtils;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.diagnostic.CliDiagnosticListener;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.diagnostic.SlimException;
import com.rtg.util.diagnostic.Talkback;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.LogFile;

import htsjdk.samtools.SAMException;
import htsjdk.samtools.util.RuntimeEOFException;
import htsjdk.samtools.util.RuntimeIOException;


/**
 * Base class for command-line modules that centralizes handling of
 * exceptions and license checking. This class does not supply any
 * logging support, for that you should extend from
 * <code>LoggedCli</code>.
 *
 */
public abstract class AbstractCli {

  private static final String LOG_DEST = System.getProperty("com.rtg.launcher.logdest", "");

  /** Usage metrics */
  protected UsageMetric mUsageMetric = null;

  /** For testing. Set each time a module is called. Reset if ever accessed. */
  private static UsageTracking sLastUsageTracking = null;
  private static void setLastUsageTracking(final UsageTracking usageTracking) {
    sLastUsageTracking = usageTracking;
  }
  /**
   * @return the log from the most recent usage tracker used.
   */
  public static String lastUsageLog() {
    if (sLastUsageTracking == null) {
      return null;
    }
    final String ret = sLastUsageTracking.usageLog();
    setLastUsageTracking(null);
    return ret;
  }

  /** The tracker for the usage metrics */
  protected UsageTracking mUsageTracker = null;

  /** Contains command configuration */
  protected CFlags mFlags;

  protected CliDiagnosticListener mMainListener;

  protected boolean mSuppressUsage;

  /**
   * @return a String with a human readable version of all the usage messages so far.
   */
  public String usageLog() {
    return mUsageTracker.usageLog();
  }

  /**
   * Initialize command line flags processor.
   * That is register the flags descriptions etc against <code>mFlags</code>
   */
  protected abstract void initFlags();

  /**
   * Perform the work of the module. Do not capture Exceptions within
   * this method -- they should be allowed to propagate up past this
   * method so they can be handled by <code>AbstractCli</code>.
   *
   * @param out standard output.
   * @param err standard error
   * @return return code.
   * @exception IOException if an I/O error occurs.
   */
  protected abstract int mainExec(final OutputStream out, final PrintStream err) throws IOException;

  /**
   * Get name of module to be used in external communications such as
   * warning messages, command line information etc.
   * @return name of application.
   */
  public abstract String moduleName();

  /**
   * Get name of application to be used in external communications such as
   * warning messages, command line information etc.
   * @return name of application.
   */
  public String applicationName() {
    return Constants.APPLICATION_NAME;
  }

  /**
   * Main program. Use -h to get help.  DO NOT call this for
   * testing. Calls <code>mainInit</code> to do the work of the
   * module, and then calls <code>System.exit</code> with the desired
   * return code.
   *
   * @param args command line arguments.
   */
  protected void mainExit(final String[] args) {
    if (!License.checkLicense()) {
      throw new NoTalkbackSlimException(ErrorType.INVALID_LICENSE);
    }
    try {
      System.exit(mainInit(args, FileUtils.getStdoutAsOutputStream(), System.err));
    } catch (final SlimException e) {
      System.exit(1);
    }
  }

  protected void createFlags(Appendable out, Appendable err) {
    final String appName = applicationName() + (moduleName().length() > 0 ? (" " + moduleName()) : "");
    mFlags = new CFlags(appName, out, err);
  }

  /**
   * Entry point for main command-line work by parsing the
   * command-line arguments, setting up listeners, calling
   * <code>mainExec</code> (implemented by subclasses), and handling
   * all Exceptions.
   *
   * @param args command line arguments
   * @param out default output (not necessarily used).
   * @param err where error and warning messages are written for user consumption.
   * @return exit code - 0 if all ok - 1 if command line arguments failed.
   */
  public int mainInit(final String[] args, final OutputStream out, final PrintStream err) {

    // Add the ability to redirect logs to a specific location, handy for non-LoggedCli commands
    switch (LOG_DEST) {
      case "":
        Diagnostic.setLogStream();
        break;
      case "err":
        Diagnostic.setLogStream(System.err);
        break;
      default:
        Diagnostic.setLogStream(new LogFile(new File(LOG_DEST)));
        break;
    }

    mUsageMetric = new UsageMetric();
    final PrintStream outPs = new PrintStream(out);
    try {
      mMainListener = initializeMainListener(err, outPs);
      try {
        final UsageTracking usage = License.usageTracking(moduleName(), CommandLine.getRunId(), mSuppressUsage);
        mUsageTracker = usage;
        setLastUsageTracking(usage);
        if (handleFlags(args, outPs, err)) {
          usage.recordBeginning();
          try {
            final int code = mainExec(out, err);
            usage.recordEnd(mUsageMetric.getMetric(), code == 0);
            if (code == 0 && !GlobalFlags.finalAccessCheck()) {
              return 1;
            }
            return code;
          } catch (final Throwable t) {
            usage.recordEnd(mUsageMetric.getMetric(), false);
            throw t;
          }
        } else {
          //mFlags.error(mFlags.getInvalidFlagMsg());
          return 1;
        }
      } catch (final IOException e) {
        Diagnostic.errorNoLog(ErrorType.IO_ERROR, getChainedErrorMessage(e));
        return 1;
      } catch (final SAMException e) {
        if (e instanceof RuntimeIOException
            || e instanceof RuntimeEOFException) {
          Diagnostic.errorNoLog(ErrorType.IO_ERROR, getChainedErrorMessage(e));
          return 1;
        }
        Diagnostic.errorNoLog(ErrorType.SAM_BAD_FORMAT_NO_FILE, e.getMessage());
        return 1;
      } catch (final SlimException e) {
        e.printErrorNoLog();
        e.invokeTalkback();
        return 1;
      } catch (final OutOfMemoryError e) {
        Diagnostic.oomMessage();
        //e.printStackTrace();
        return 1;
      } catch (final Throwable t) { //catch everything except SlimException which has already been logged etc.
        err.println(t.getMessage());
        Diagnostic.errorNoLog(ErrorType.SLIM_ERROR);
        if (License.isDeveloper()) {
          t.printStackTrace(err);
        }
        Talkback.postTalkback(t);
        return 1;
      } finally {
        Diagnostic.removeListener(mMainListener);
      }
    } finally {
      outPs.flush();
    }
  }

  private String getChainedErrorMessage(Throwable e) {
    final StringBuilder message = new StringBuilder();
    message.append(e.getMessage());
    for (Throwable t : e.getSuppressed()) {
      message.append(StringUtils.LS).append(t.getMessage());
    }
    return message.toString();
  }

  /**
   * @param err errors and warnings for user, also final last resort place to put logs if cannot create one in working directory.
   * @param out information messages for user
   * @return listener set
   */
  protected CliDiagnosticListener initializeMainListener(PrintStream err, PrintStream out) {
    final CliDiagnosticListener listener = new CliDiagnosticListener(err, out);
    Diagnostic.addListener(listener);
    return listener;
  }

  /**
   * Create the flags object.
   * @param out output stream
   * @param err error stream
   */
  protected void createRegisterFlags(Appendable out, Appendable err) {
    createFlags(out, err);
    initFlags();
  }

  //separate so can use in tests
  protected boolean handleFlags(String[] args, Appendable out, Appendable err) {
    createRegisterFlags(out, err);
    if (!GlobalFlags.initialAccessCheck()) {
      return false;
    }
    GlobalFlags.registerExperimentalFlags(mFlags);
    return mFlags.setFlags(args);
  }

  /**
   * Get the flags.
   * @return command line flags
   */
  protected CFlags getCFlags() {
    return mFlags;
  }

}
