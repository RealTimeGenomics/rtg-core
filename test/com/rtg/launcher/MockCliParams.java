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

import static com.rtg.launcher.BuildCommon.MAXGAP_FLAG;

import java.io.File;

import com.rtg.util.InvalidParamsException;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.Validator;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;


/**
 * Used for testing error handling etc. in <code>CLI</code>.
 * Will throw the errors depending on the flags set.
 */
public class MockCliParams extends ModuleParams {


  private static final class MockValidator implements Validator {
    @Override
    public boolean isValid(final CFlags flags) {
      if (flags.isSet("validator")) {
        Diagnostic.error(ErrorType.INVALID_INTEGER_FLAG_VALUE, MAXGAP_FLAG, Integer.toString(42));
        return false;
      }
      return true;
    }
  }

  /**
   * Set flags.
   * @param flags to be set.
   */
  public static void makeFlags(final CFlags flags) {
    flags.setValidator(new MockValidator());
    flags.registerOptional('g', "global", "error thrown in globalIntegrity()");
    flags.registerOptional('v', "validator", "error thrown in validator");
    flags.registerOptional('c', "constructor", "InvalidParamsException thrown in params constructor");
    //this happens before the log is switched - it should be left there and include the cause
    flags.registerOptional('x', "constructorx", "RuntimeException thrown in params constructor");
    flags.registerOptional('e', "runtimeErr", "Write to err during task execution");
    flags.registerOptional('r', "runtime", "RuntimeException thrown during execution of task");
    flags.registerOptional('s', "runtimeslim", "SlimException thrown during execution of task");
  }

  private final boolean mGlobalError;

  private final boolean mValidatorError;

  private final boolean mConstructorError;

  final boolean mRuntime;

  final boolean mRuntimeSlim;

  final boolean mRuntimeErr;

  /**
   * @param flags command line flags.
   * @throws InvalidParamsException if there are errors in the command line parameters.
   */
  public MockCliParams(final CFlags flags) {
    super(flags.getName());
    mGlobalError = flags.isSet("global");
    mValidatorError = flags.isSet("validator");
    mConstructorError = flags.isSet("constructor");
    final boolean errorX = flags.isSet("constructorx");
    mRuntime = flags.isSet("runtime");
    mRuntimeSlim = flags.isSet("runtimeslim");
    mRuntimeErr = flags.isSet("runtimeErr");

    if (mConstructorError) {
      throw new InvalidParamsException(ErrorType.INVALID_MAX_INTEGER_FLAG_VALUE, "-c", "42", "41");
    }
    if (errorX) {
      throw new RuntimeException("in parameter constructor");
    }

  }

  @Override
  public String toString() {
    return (mGlobalError ? " global" : "") + (mValidatorError ? " validator" : "") + (mConstructorError ? " constructor" : "");
  }

  @Override
  public File file(final String name) {
    throw new UnsupportedOperationException();
  }

  private File mDirectory = null;

  public void setDirectory(final File dir) {
    mDirectory = dir;
  }

  @Override
  public File directory() {
    if (mDirectory != null) {
      return mDirectory;
    }
    throw new UnsupportedOperationException();
  }

}
