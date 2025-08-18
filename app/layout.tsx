import "./globals.css";
import type { Metadata } from "next";
import Navbar from "../components/Navbar";
import ProgressProvider from "../components/ProgressProvider";
import ProgressPane from "../components/ProgressPane";

export const metadata: Metadata = {
  title: "EMTR UI",
  description: "Run EMTR and view results",
};

export default function RootLayout({ children }: { children: React.ReactNode }) {
  return (
    <html lang="en">
      <body className="bg-gray-500 text-white">
        <Navbar />
        <ProgressProvider>
          <div className="max-w-6xl mx-auto px-4">
            <div className="flex flex-col lg:flex-row gap-4 py-4">
              <main className="flex-1">{children}</main>
              <ProgressPane />
            </div>
          </div>
        </ProgressProvider>
      </body>
    </html>
  );
}
