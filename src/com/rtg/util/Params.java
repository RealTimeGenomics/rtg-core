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
package com.rtg.util;

import java.io.Closeable;

/**
 * An interface implemented by most of the parameter classes.
 * The <a href="http://twit88.com/blog/2008/08/31/java-create-immutable-object-using-builder-pattern/">Builder Pattern</a> should be used for their construction.
 */
public interface Params extends Closeable {

  /**
   * Get status of underlying files.
   * @return true iff all underlying file objects are closed.
   */
  boolean closed();

}
