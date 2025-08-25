"use client";

import React from "react";

export default function SignificanceOfRecallPage() {
  // cache-bust param so a fresh image shows after you regenerate the PNG
  const imgSrc = "/analysis/recall_distribution.PNG?v=" + Date.now();

  return (
    <main className="px-6 py-8 max-w-5xl mx-auto">
      <h1 className="text-3xl font-semibold mb-2">Significance of Recall</h1>
      <p className="text-sm text-gray-600 mb-6">
        Distribution of recall values obtained by rooting the unrooted tree at internal nodes.
      </p>

      <div className="rounded-2xl shadow border p-4 bg-white">
        {/* Plain <img> avoids needing fixed width/height for Next/Image */}
        <img
          src={imgSrc}
          alt="Distribution of recall values"
          className="w-full h-auto rounded-lg"
        />
      </div>
    </main>
  );
}