// app/violin/page.tsx
'use client';
import React from "react";

export default function ViolinPage() {
  return (
    <div className="p-6">
      <h1 className="text-2xl font-bold mb-4">Violin Plots</h1>
      <p className="text-sm text-gray-600">
        Select a dataset and render violin plots here.
      </p>
      {/* TODO: mount your violin plot component here */}
    </div>
  );
}
