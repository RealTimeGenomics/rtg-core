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

import java.util.Collection;

/**
 * Create Collections for backing the multi-map.
 * @param <W> type of the collection.
 */
public interface MultiMapFactory <W> {
    /**
     * Create a new collection used for each instance of the key.
     * @return a new collection.
     */
    Collection<W> createCollection();
}
