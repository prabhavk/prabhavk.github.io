// app/layout.tsx
import "./globals.css";
import type { Metadata } from "next";
import Nav from "@/components/Nav"; 

export const metadata: Metadata = {
  title: "emtree UI",
  description: "Run emtree and view results",
};

export default function RootLayout({ children }: { children: React.ReactNode }) {
  return (
    <html lang="en">
      <body className="min-h-screen bg-gray-50 text-gray-900">
        <a
          href="#main"
          className="sr-only focus:not-sr-only focus:absolute focus:top-2 focus:left-2 bg-white border rounded px-3 py-2"
        >
          Skip to content
        </a>

        <Nav />

        <div className="max-w-7xl mx-auto pl-4 pr-0">
          <main id="main" className="py-4">{children}</main>
        </div>
      </body>
    </html>
  );
}