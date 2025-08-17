"use client";

export default function FlowerPage() {
  return (
    <main>
      <h1 className="text-2xl font-bold mb-4">Flower Plots</h1>
      <p className="text-gray-700">
        This is a placeholder page. Later, render petal/rose charts comparing methods per replicate.
      </p>
      {/* TODO:
          - Load per-replicate deltas
          - Render flower/rose plots (e.g., Plotly barpolar) */}
    </main>
  );
}
