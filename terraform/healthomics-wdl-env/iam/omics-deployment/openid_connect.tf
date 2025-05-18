resource "aws_iam_openid_connect_provider" "github" {
  count = var.aws_region == "us-east-1" ? 1 : 0 # create only once on primary region
  client_id_list = [
    "https://github.com/${var.github_org}",
    "sts.amazonaws.com"
  ]

  thumbprint_list = ["6938fd4d98bab03faadb97b34396831e3780aea1"]
  url = "https://token.actions.githubusercontent.com"
}