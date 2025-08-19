import { connect } from "@planetscale/database";

export const db = () =>
  connect({
    host: process.env.PS_HOST!,
    username: process.env.PS_USERNAME!,
    password: process.env.PS_PASSWORD!,
  });
