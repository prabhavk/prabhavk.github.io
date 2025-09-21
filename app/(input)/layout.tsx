// app/(input)/layout.tsx
"use client";

import React from "react";
import ProgressProvider from "@/components/ProgressProvider";
import ClockAndTree from "@/components/ClockAndTree";
import DebugProgressPane from "@/components/DebugProgressPane";

export default function InputLayout({ children }: { children: React.ReactNode }) {
  
  const [showDebug, setShowDebug] = React.useState<boolean>(false);

  React.useEffect(() => {
    try {
      const saved = localStorage.getItem("emtr.showDebugPane");
      if (saved === "1") setShowDebug(true);
      if (saved === "0") setShowDebug(false);
    } catch {}
  }, []);

  const handleToggle = (v: boolean) => {
    setShowDebug(v);
    try {
      localStorage.setItem("emtr.showDebugPane", v ? "1" : "0");
    } catch {}
  };

  return (
    <ProgressProvider>
      <div className="relative left-1/2 right-1/2 -ml-[50vw] -mr-[50vw] w-screen">
        <div
          className="
            grid grid-cols-1
            lg:grid-cols-[minmax(0,2fr)_minmax(0,3fr)]
            xl:grid-cols-[minmax(0,1.6fr)_minmax(0,2.4fr)]
            gap-3
            pr-0 pl-6 sm:pl-8 lg:pl-10 xl:pl-12 2xl:pl-16
          "
        >
          {/* Input form goes here */}
          <main id="main" className="min-w-0">
            {children}
          </main>

          {/* Sidebar: toggle between DebugProgressPane and ClockAndTree */}
          <aside className="lg:sticky lg:top-4 self-start max-h-[calc(100vh-0.5rem)] overflow-auto">
            <div className="flex items-center justify-between mb-2 px-1">              
              <label className="flex items-center gap-2 text-xs text-black">
                <span>Show debug logs</span>
                <input
                  type="checkbox"
                  checked={showDebug}
                  onChange={(e) => handleToggle(e.target.checked)}
                  className="accent-black"
                  aria-label="Toggle debug progress pane"
                />
              </label>
            </div>

            {showDebug ? (
              <DebugProgressPane height={550} />
            ) : (
              <ClockAndTree clockHeight={400} treeHeight={500} yShift={-50}/>
            )}
          </aside>
        </div>
      </div>
    </ProgressProvider>
  );
}
