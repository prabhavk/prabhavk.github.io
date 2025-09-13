// app/cake/page.tsx
"use client";

import React, { useEffect, useState } from "react";

export default function TreesPage() {
  const [job, setJob] = useState<string>("");

  // Load job from URL or localStorage
  useEffect(() => {
    try {
      const u = new URL(window.location.href);
      const fromUrl = u.searchParams.get("job") || u.searchParams.get("job_id") || "";
      const fromLs = localStorage.getItem("emtr:selectedJobId") || "";
      setJob(fromUrl || fromLs || "");
    } catch {}
  }, []);

  const imageUrl = "/assets/evergreen.png"; // strongly watermarked version

  return (
    <div className="relative min-h-screen text-white overflow-hidden">
      {/* Background image with filters to hide watermark */}
      <div
        className="absolute inset-0"
        style={{
          backgroundImage: `url(${imageUrl})`,
          backgroundSize: "cover",
          backgroundPosition: "center",
          backgroundRepeat: "no-repeat",
          filter: "contrast(90%) brightness(115%) blur(0.3px)", // soften/hide watermark
          opacity: 0.75, // fade overall image slightly
        }}
      />

      {/* Optional subtle white veil to hide watermark further */}
      <div className="absolute inset-0 bg-white/15" />

      {/* Foreground content */}
      <div className="relative z-10 max-w-6xl mx-auto p-6">
        <div className="text-sm text-gray-200">
          Job: <span className="font-mono">{job || "(none)"}</span>
        </div>
      </div>
    </div>
  );
}
