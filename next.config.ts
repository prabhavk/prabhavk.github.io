import type { NextConfig } from "next";

const devExperimental =
  process.env.NODE_ENV === "development"
    ? {
        experimental: {
          // allow your LAN device to open the dev server
          allowedDevOrigins: ["http://192.168.29.107:3000"],
        },
      }
    : {};

const nextConfig: NextConfig = {
  output: "export",
  ...devExperimental,
};

export default nextConfig;
